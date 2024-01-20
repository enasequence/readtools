package uk.ac.ebi.ena.readtools.v2;

import uk.ac.ebi.ena.readtools.v2.provider.ReadsProvider;
import uk.ac.ebi.ena.readtools.v2.provider.ReadsProviderFactory;
import uk.ac.ebi.ena.readtools.v2.read.IRead;
import uk.ac.ebi.ena.readtools.v2.validator.ReadsValidationException;

public class MockReadsProviderFactory extends ReadsProviderFactory {
    private final ReadsProvider<? extends IRead> provider;

    public MockReadsProviderFactory(MockReadsProvider.MockRead... reads) {
        super(null, null);
        this.provider = new MockReadsProvider(reads);
    }

    @Override
    public ReadsProvider<? extends IRead> makeReadsProvider() throws ReadsValidationException {
        return provider;
    }
}
