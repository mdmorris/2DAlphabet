---
layout: default
title: Configuration files
has_children: true
nav_order: 3
---

# Configuration files

The goal of the JSON is to have an easily understandable and configurable input
that allows the
2D Alphabet software to read and organize analysis files and histograms to its
liking while also giving the user the ability to easily configure values like
bin sizes and ranges without lots of command line options. This means that while
the user is encouraged to add what they need, there are some keys and values
that must stay the same. These static strings are always in capital letters
to make them easy to distinguish. The six static sections are described below.
