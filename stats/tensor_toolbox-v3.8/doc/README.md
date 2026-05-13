# Documentation Instructions

In general, HTML documentation is generated from files of the form
`doc/XXX_doc.m` where `XXX` is the name of the function or class.

To create the HTML documentation, be sure to use the style sheet 
`ttb.xsl` as follows:

```
cd doc
publish('XXX_doc.m','stylesheet','ttb.xsl') 
```

This will create the file `doc/html/XXX_doc.html` and any accompanying 
images.

## Important

- The code should be 100% reproducible and results should be exactly the
  same each time.  This means that any random number generation should be 
  initialized with `rng('default')` or `rng(seed)` at the beginning of 
  each code cell.
- Give each problem a unique name so a single part of the code can be
  re-run and still produce the same results.

# Special Files

There are a few special files that may also need edits.

1. `doc/html/helptoc.xml` is the table of contents for the HTML 
	documentation which appears on the left-hand side. 
2. HTML Files not generated from `doc/XXX_doc.m`:
   - `doc/html/index.html` is the main page for the Tensor Toolbox.
   - `doc/html/tensor_types.html` is for tensor types. 
   - `doc/html/cp.html` is for CP decomposition methods. 
   - `doc/html/tucker.html` is for Tucker decomposition methods. 
   - `doc/html/eigen.html` is for tensor eigenvalue methods.
   - `doc/html/working_with_tensors.html` is for working with tensors.
   - `doc/html/converting.html` links to the tensor as matrices classes.

# Notes

When the documentation is viewed via MATLAB's `doc` command, the page created is viewed inside an iframe, and the contents links don't exactly work properly.  To see the documentation as it would appear on the web, open the file `doc/html/index.html` in a web browser.