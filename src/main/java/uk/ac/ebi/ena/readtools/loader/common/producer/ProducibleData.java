/*
* Copyright 2010-2020 EMBL - European Bioinformatics Institute
* Licensed under the Apache License, Version 2.0 (the "License"); you may not use this
* file except in compliance with the License. You may obtain a copy of the License at
* http://www.apache.org/licenses/LICENSE-2.0
* Unless required by applicable law or agreed to in writing, software distributed under the
* License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
* CONDITIONS OF ANY KIND, either express or implied. See the License for the
* specific language governing permissions and limitations under the License.
*/
package uk.ac.ebi.ena.readtools.loader.common.producer;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * A member field of a class annotated with this is expected to hold data that will be acquired by calling the method
 * whose name is defined in the property {@link ProducibleData#method}.
 */
@Retention( RetentionPolicy.RUNTIME )
@Target( {ElementType.FIELD } )
public @interface
ProducibleData
{
    /**
     * The name of the method that is supposed to return the data to be held by this field.
     */
    String method();
}