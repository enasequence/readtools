/*
* Copyright 2010-2021 EMBL - European Bioinformatics Institute
* Licensed under the Apache License, Version 2.0 (the "License"); you may not use this
* file except in compliance with the License. You may obtain a copy of the License at
* http://www.apache.org/licenses/LICENSE-2.0
* Unless required by applicable law or agreed to in writing, software distributed under the
* License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
* CONDITIONS OF ANY KIND, either express or implied. See the License for the
* specific language governing permissions and limitations under the License.
*/
package uk.ac.ebi.ena.readtools.loader.common;

import java.util.Collection;

public class InvalidBaseCharacterException extends RuntimeException {
    private final String bases;

    private final Collection<Character> invalidCharacters;

    public InvalidBaseCharacterException(String bases, Collection<Character> invalidCharacters) {
        this.bases = bases;
        this.invalidCharacters = invalidCharacters;
    }

    public String getBases() {
        return bases;
    }

    public Collection<Character> getInvalidCharacters() {
        return invalidCharacters;
    }

    @Override
    public String getMessage() {
        return String.format("Invalid characters %s in bases : %s", invalidCharacters, bases);
    }
}
