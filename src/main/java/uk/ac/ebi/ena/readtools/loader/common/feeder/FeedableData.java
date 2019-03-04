package uk.ac.ebi.ena.readtools.loader.common.feeder;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

@Retention( RetentionPolicy.RUNTIME )
@Target( {ElementType.FIELD } )
public @interface
FeedableData
{
    String method();
}