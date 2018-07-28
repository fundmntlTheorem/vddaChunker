# vddaChunker
Tools for evaluating our chunking algorithms
These are a set of files that we can use to make an even playing field for testing out our ideas for chunking a DDA spectrum.  Require the module chunkingEvaluator.lua, and LoadSpectra to get a set of spectra with peaks to use.  For each spectrum, create a list of chunks, which are tables with {mass = center mass, width = isolation full width}.  Use EvaluateAChunking( spectrum, Chunker ) to evaluate one spectrum, or EvaluateFullChunkingJob( spectra, Chunker ) to evaluate your job on a set of spectra.
