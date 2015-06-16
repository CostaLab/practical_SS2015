# What

First, all the motifs of a specific technology, for all organisms, were *joined together*.

Then, common motifs were calculated for different technology combinations.

# Preconditions

Since different organisms' results are joined together (by technology), some results might be uncorrect.

For example, immagine a situation where, using GAIIx, we have similar motifs for a few organisms but
a widely different set of motifs for, say, Bordetella Pertussis.

Also, suppose that GAIIx and MiSeq normally vary a lot in motifs found.

By joining all organisms together for GAIIx, we are much more likely to detect common motifs between
GAIIx and MiSeq even though they are due to an outlier (B. Pertussis)!

This is just a fictional example, but the possibility is there. So any interpretation of data should
also rely on analysis of similarity for a specific sequencer among the different organisms.

# Why

Our assumption is that motifs for different platforms are not necessarily consistent - possibly not at all.

Some technologies might yeld similar motifs though, ie Illumina's sequencers have probably similar motifs.

# Results

TODO
