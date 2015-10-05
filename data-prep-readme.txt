Documentation on how to get lprobs object from data.

Writing this out here to record the manual process.
First get the data to a clean table form.

Then, using limma:
.
.
.
Single-vs-WT log-probs should be recorded in a reporters x actors matrix.
Double-vs-Single log-probs should be recorded in a list (indexed by actors) of lists (indexed by actors) of vectors (of length reporters).