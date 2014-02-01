NOTES ON MEX INTERFACE
----------------------


(1) Pitfalls:

Quick and dirty interface!

- Macros M_GET*, M_SET*, M_MAKE* depend on argidx (int). Each starts with
  ++argidx, so set
    argidx = -1;
  to start from the beginning. Has to be done separately for M_GET* vs.
  M_SET*, M_MAKE*

- Array variables:
    T* xxx; int nxxx; // T is int or double
  Pass to wrapper: M_ARR(xxx)
  NOTE: We do not in general use fst_vector.

- Create array:
  When calling M_MAKE?ARRAY(xxx), nxxx must contain the size (and argidx be
  set appropriately)

- Macros for exported functions:
  Most exported functions have 'char* errstr' as last argument, has to point
  to error message buffer (enough space: 512). These are defined as macros
  XXX(...), the real functions are _XXX(...,char* errstr)
