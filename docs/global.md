---
layout: default
title: GLOBAL
has_children: false
parent: Configuration files
nav_order: 5
---

# `GLOBAL`
This section is designed to help users with large configuration file
names by allowing them to create JSON-wide variables. For example,
if all of your files are located in `/long/path/to/my/files/`, you 
can store this string in the GLOBAL dictionary with a custom key 
(let's say `path`). Now instead of having to write the full directory
path for every process and systematic, the user can just write `path`.
This simplifies the JSON and also has the standard advantages of using
variables over several instances of the same object.

This functionality works by searching all strings in the JSON for instances
of each key in `GLOBAL` and replacing the key with its corresponding dictionary value.
Note that it does not work for booleans or integers/floats.

The user must be careful they don't accidentally use strings in the JSON
that are identical to keys in `GLOBAL` so accidental substitutions don't happen.
This means keys in `GLOBAL` should be at least partially descriptive 
(single character keys would be a bad idea). 