<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">

<html xmlns="http://www.w3.org/1999/xhtml" lang="en" xml:lang="en">
<head>
    <meta http-equiv="Content-Type" content="text/html; charset=us-ascii" />
  <title>Boost C++ Libraries - boost/utility/swap.hpp</title>  <link rel="icon" href="/favicon.ico" type="image/ico" />
  <link rel="stylesheet" type="text/css" href="/style/section-doc.css" />
  <!--[if IE]> <style type="text/css"> body { behavior: url(/style/csshover.htc); } </style> <![endif]-->

</head><!-- boost_1_39_0/boost/utility/swap.hpp -->

<body>
  <div id="heading">
      <div id="heading-placard"></div>

  <h1 id="heading-title"><a href="/"><img src="/gfx/space.png" alt=
  "Boost C++ Libraries" id="heading-logo" /><span id="boost">Boost</span>
  <span id="cpplibraries">C++ Libraries</span></a></h1>

  <p id="heading-quote"><span class="quote">&ldquo;...one of the most highly
  regarded and expertly designed C++ library projects in the
  world.&rdquo;</span> <span class="attribution">&mdash; <a href=
  "http://www.gotw.ca/" class="external">Herb Sutter</a> and <a href=
  "http://en.wikipedia.org/wiki/Andrei_Alexandrescu" class="external">Andrei
  Alexandrescu</a>, <a href=
  "http://safari.awprofessional.com/?XmlId=0321113586" class="external">C++
  Coding Standards</a></span></p>

  <div id="heading-sections">
    <ul>
      <li id="welcome-section-tab"><a href="/">Welcome</a></li>

      <li id="boost-section-tab"><a href="/users/">Introduction</a></li>

      <li id="community-section-tab"><a href="/community/">Community</a></li>

      <li id="development-section-tab"><a href=
      "/development/">Development</a></li>

      <li id="support-section-tab"><a href="/support/">Support</a></li>

      <li id="doc-section-tab"><a href="/doc/">Documentation</a></li>

      <li id="map-section-tab"><a href="/map.html">Index</a></li>
    </ul>
  </div>
  </div>

  <div id="body">
    <div id="body-inner">
      <div id="content">
        <div class="section" id="docs">
          <div class="section-0">
            <div class="section-body">
              <h3>boost/utility/swap.hpp</h3>
<pre>
// Copyright (C) 2007, 2008 Steven Watanabe, Joseph Gauterin, Niels Dekker
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
// For more information, see http://www.boost.org


#ifndef BOOST_UTILITY_SWAP_HPP
#define BOOST_UTILITY_SWAP_HPP

// Note: the implementation of this utility contains various workarounds:
// - swap_impl is put outside the boost namespace, to avoid infinite
// recursion (causing stack overflow) when swapping objects of a primitive
// type.
// - swap_impl has a using-directive, rather than a using-declaration,
// because some compilers (including MSVC 7.1, Borland 5.9.3, and
// Intel 8.1) don't do argument-dependent lookup when it has a
// using-declaration instead.
// - boost::swap has two template arguments, instead of one, to
// avoid ambiguity when swapping objects of a Boost type that does
// not have its own boost::swap overload.

#include &lt;algorithm&gt; //for std::swap
#include &lt;cstddef&gt; //for std::size_t

namespace boost_swap_impl
{
  template&lt;class T&gt;
  void swap_impl(T&amp; left, T&amp; right)
  {
    using namespace std;//use std::swap if argument dependent lookup fails
    swap(left,right);
  }

  template&lt;class T, std::size_t N&gt;
  void swap_impl(T (&amp; left)[N], T (&amp; right)[N])
  {
    for (std::size_t i = 0; i &lt; N; ++i)
    {
      ::boost_swap_impl::swap_impl(left[i], right[i]);
    }
  }
}

namespace boost
{
  template&lt;class T1, class T2&gt;
  void swap(T1&amp; left, T2&amp; right)
  {
    ::boost_swap_impl::swap_impl(left, right);
  }
}

#endif
</pre>
            </div>
          </div>
        </div>
      </div>

      <div class="clear"></div>
    </div>
  </div>

  <div id="footer">
    <div id="footer-left">
      <div id="revised">
        <p>Revised $Date: 2008-11-03 08:44:40 -0500 (Mon, 03 Nov 2008) $</p>
      </div>

      <div id="copyright">
        <p>Copyright Beman Dawes, David Abrahams, 1998-2005.</p>

        <p>Copyright Rene Rivera 2004-2008.</p>
      </div>  <div id="license">
    <p>Distributed under the <a href="/LICENSE_1_0.txt" class=
    "internal">Boost Software License, Version 1.0</a>.</p>
  </div>
    </div>

    <div id="footer-right">
        <div id="banners">
    <p id="banner-xhtml"><a href="http://validator.w3.org/check?uri=referer"
    class="external">XHTML 1.0</a></p>

    <p id="banner-css"><a href=
    "http://jigsaw.w3.org/css-validator/check/referer" class=
    "external">CSS</a></p>

    <p id="banner-osi"><a href=
    "http://www.opensource.org/docs/definition.php" class="external">OSI
    Certified</a></p>
  </div>
    </div>

    <div class="clear"></div>
  </div>
</body>
</html>
