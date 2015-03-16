## Format ##

Put your own normalization constants in a text file such as "norm.txt" and provide it to diffReps using the "--norm" option. The program will only read the first two lines of the normalization file and ignore the rest. Make sure you separate the group name and normalization constants by white space. Each group must be on a different line. An example is as follows:

```
treatment       0.9 1.1 1.2
control         1.3 1.5 0.6
```

everything else is ignored here...