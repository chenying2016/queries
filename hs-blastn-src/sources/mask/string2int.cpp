#include "def.h"
#include "mask_macros.h"

typedef string CTempString;

    /// String to number conversion flags.
    enum EStringToNumFlags {
        fConvErr_NoThrow      = (1 << 16),   ///< On error, return zero and set
        /// errno to non-zero instead of throwing an exception (the default).
        /// We recommend the following technique to check against errors
        /// with minimum overhead when this flag is used:
        ///     if (!retval  &&  errno != 0)
        ///        ERROR;
        /// And for StringToDouble*() variants:
        ///     if (retval == HUGE_VAL  ||  retval == -HUGE_VAL  ||
        ///        !retval  &&  errno != 0)
        ///        ERROR;

        fMandatorySign        = (1 << 17),   ///< See 'fWithSign'
        fAllowCommas          = (1 << 18),   ///< See 'fWithCommas'
        fAllowLeadingSpaces   = (1 << 19),   ///< Can have leading spaces
        fAllowLeadingSymbols  = (1 << 20) | fAllowLeadingSpaces,
                                             ///< Can have leading non-nums
        fAllowTrailingSpaces  = (1 << 21),   ///< Can have trailing spaces
        fAllowTrailingSymbols = (1 << 22) | fAllowTrailingSpaces,
                                             ///< Can have trailing non-nums
        fDecimalPosix         = (1 << 23),   ///< For decimal point, use C locale
        fDecimalPosixOrLocal  = (1 << 24),   ///< For decimal point, try both C and current locale

        fDS_ForceBinary       = (1 << 26),
        fDS_ProhibitFractions = (1 << 27),
        fDS_ProhibitSpaceBeforeSuffix = (1 << 28)
    };
    typedef int TStringToNumFlags;   ///< Bitwise OR of "EStringToNumFlags"

// Skip all allowed chars (all except used for digit composition).
// Update 'ptr' to current position in the string.
enum ESkipMode {
    eSkipAll,           // all symbols
    eSkipAllAllowed,    // all symbols, except digit/+/-/.
    eSkipSpacesOnly     // spaces only
};

static inline
bool s_IsDecimalPoint(unsigned char ch, TStringToNumFlags  flags)
{
    if ( ch != '.' && ch != ',') {
        return false;
    }
    if (flags & fDecimalPosix) {
        return ch == '.';
    }
    else if (flags & fDecimalPosixOrLocal) {
        return ch == '.' || ch == ',';
    }
    struct lconv* conv = localeconv();
    return ch == *(conv->decimal_point);
}

#define CHECK_COMMAS                                                  \
/* Check on possible commas */                                    \
    if (flags & fAllowCommas) {                                 \
        if (ch == ',') {                                              \
		if ((numpos == pos)  ||                                   \
			((comma >= 0)  &&  (comma != 3)) ) {                  \
			/* Not first comma, sitting on incorrect place */     \
			break;                                                \
		}                                                         \
		/* Skip it */                                             \
		comma = 0;                                                \
		pos++;                                                    \
		continue;                                                 \
        } else {                                                      \
            if (comma >= 0) {                                         \
			/* Count symbols between commas */                    \
			comma++;                                              \
		}                                                         \
	}                                                             \
}

/// @internal
// Check that symbol 'ch' is good symbol for number with radix 'base'.
static inline
bool s_IsGoodCharForRadix(char ch, int base, int* value = 0)
{
    if ( base <= 10 ) {
        // shortcut for most frequent case
        int delta = ch-'0';
        if ( unsigned(delta) < unsigned(base) ) {
            if ( value ) {
                *value = delta;
            }
            return true;
        }
        return false;
    }
    if (!isalnum((unsigned char) ch)) {
        return false;
    }
    // Corresponding numeric value of *endptr
    int delta;
    if (isdigit((unsigned char) ch)) {
        delta = ch - '0';
    } else {
        ch = tolower((unsigned char) ch);
        delta = ch - 'a' + 10;
    }
    if ( value ) {
        *value = delta;
    }
    return delta < base;
 }

static inline
void s_SkipAllowedSymbols(const CTempString& str,
                          SIZE_TYPE&         pos,
                          ESkipMode          skip_mode,
                          TStringToNumFlags  flags)
{
    if (skip_mode == eSkipAll) {
        pos = str.length();
        return;
    }

    for ( SIZE_TYPE len = str.length(); pos < len; ++pos ) {
        unsigned char ch = str[pos];
        if ( isdigit(ch)  ||  ch == '+' ||  ch == '-'  ||  s_IsDecimalPoint(ch,flags) ) {
            break;
        }
        if ( (skip_mode == eSkipSpacesOnly)  &&  !isspace(ch) ) {
            break;
        }
    }
}

static inline
bool s_CheckRadix(const CTempString& str, SIZE_TYPE& pos, int& base)
{
    if ( base == 10 || base == 8 ) {
        // shortcut for most frequent case
        return true;
    }
    // Check base
    if ( base < 0  ||  base == 1  ||  base > 36 ) {
        return false;
    }
    // Try to determine base using first chars of the string
    unsigned char ch   = str[pos];
    unsigned char next = str[pos+1];
    if ( base == 0 ) {
        if ( ch != '0' ) {
            base = 10;
        } else if (next == 'x' || next == 'X') {
            base = 16;
        } else {
            base = 8;
        }
    }
    // Remove leading '0x' for hex numbers
    if ( base == 16 ) {
        if (ch == '0'  &&  (next == 'x' || next == 'X')) {
            pos += 2;
        }
    }
    return true;
}


Uint8 StringToUInt8(const CTempString& str,
                    TStringToNumFlags flags, int base)
{
    ASSERT(flags == 0  ||  flags > 32);

    const TStringToNumFlags slow_flags =
        fMandatorySign|fAllowCommas|fAllowLeadingSymbols|fAllowTrailingSymbols;

    if ( base == 10 && (flags & slow_flags) == 0 ) {
        // fast conversion

        // Current position in the string
        CTempString::const_iterator ptr = str.begin(), end = str.end();

        // Determine sign
        if ( ptr != end && *ptr == '+' ) {
            ++ptr;
        }
        //if ( ptr == end ) {
        //    S2N_CONVERT_ERROR(Uint8, kEmptyStr, EINVAL, ptr-str.begin());
        //}
        ASSERT(ptr != end);

        // Begin conversion
        Uint8 n = 0;

        const Uint8 limdiv = kMax_UI8/10;
        const int   limoff = int(kMax_UI8 % 10);

        do {
            char ch = *ptr;
            int  delta = ch - '0';
            //if ( unsigned(delta) >= 10 ) {
            //    S2N_CONVERT_ERROR(Uint8, kEmptyStr, EINVAL, ptr-str.begin());
            //}
            ASSERT(unsigned(delta) < 10);
            // Overflow check
            if ( n >= limdiv && (n > limdiv || delta > limoff) ) {
                //S2N_CONVERT_ERROR(Uint8, kEmptyStr, ERANGE, ptr-str.begin());
            }
            n = n*10+delta;
        } while ( ++ptr != end );

        return n;
    }

    // Current position in the string
    SIZE_TYPE pos = 0, size = str.size();

    // Skip allowed leading symbols
    if (flags & fAllowLeadingSymbols) {
        bool spaces = ((flags & fAllowLeadingSymbols) == fAllowLeadingSpaces);
        s_SkipAllowedSymbols(str, pos,
                             spaces ? eSkipSpacesOnly : eSkipAllAllowed, flags);
    }
    // Determine sign
    if (str[pos] == '+') {
        pos++;
    } else {
        if (flags & fMandatorySign) {
            //S2N_CONVERT_ERROR_INVAL(Uint8);
        }
    }
    SIZE_TYPE pos0 = pos;

    // Begin conversion
    Uint8 n = 0;
    // Check radix base
    if ( !s_CheckRadix(str, pos, base) ) {
        //S2N_CONVERT_ERROR_RADIX(Uint8, "bad numeric base '" +
        //                        NStr::IntToString(base) + "'");
    }

    Uint8 limdiv = kMax_UI8 / base;
    int   limoff = int(kMax_UI8 % base);

    // Number of symbols between two commas. '-1' means -- no comma yet.
    int       comma  = -1;
    SIZE_TYPE numpos = pos;

    while (char ch = str[pos]) {
        int delta;  // corresponding numeric value of 'ch'

        // Check on possible commas
        CHECK_COMMAS;
        // Sanity check
        if ( !s_IsGoodCharForRadix(ch, base, &delta) ) {
            break;
        }
        // Overflow check
        if ( n >= limdiv  &&  (n > limdiv  ||  delta > limoff) ) {
            //S2N_CONVERT_ERROR_OVERFLOW(Uint8);
        }
        n *= base;
        n += delta;
        pos++;
    }

    // Last checks
    if ( pos == pos0  || ((comma >= 0)  &&  (comma != 3)) ) {
        //S2N_CONVERT_ERROR_INVAL(Uint8);
    }
    // Skip allowed trailing symbols
    if (flags & fAllowTrailingSymbols) {
        bool spaces = ((flags & fAllowTrailingSymbols) ==
                       fAllowTrailingSpaces);
        s_SkipAllowedSymbols(str, pos, spaces ? eSkipSpacesOnly : eSkipAll, flags);
    }
    //CHECK_ENDPTR_SIZE(Uint8);
    return n;
}

unsigned int
StringToUInt(const CTempString& str, TStringToNumFlags flags, int base)
{
    //S2N_CONVERT_GUARD_EX(flags);
    Uint8 value = StringToUInt8(str, flags, base);
    //if ( value > UINT4_MAX ) {
        //S2N_CONVERT_ERROR(unsigned int, "overflow", ERANGE, 0);
    //}
	ASSERT(value <= UINT4_MAX);
    return (unsigned int) value;
}
