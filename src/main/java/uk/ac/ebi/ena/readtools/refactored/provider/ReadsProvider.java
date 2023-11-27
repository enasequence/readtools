package uk.ac.ebi.ena.readtools.refactored.provider;

import uk.ac.ebi.ena.readtools.refactored.read.IRead;

public interface ReadsProvider<T extends IRead> extends Iterable<T>, AutoCloseable {}
