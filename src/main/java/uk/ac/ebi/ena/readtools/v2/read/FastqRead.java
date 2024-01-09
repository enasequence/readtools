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
package uk.ac.ebi.ena.readtools.v2.read;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.ebi.ena.readtools.loader.common.writer.ReadWriterException;

public class FastqRead implements IRead {
    /*
    @ Each sequence identifier line starts with @
    1    <instrument> Characters
    allowed:
    a-z, A-Z, 0-9 and
    underscore
    2    Instrument ID
    <run number> Numerical Run number on instrument
    3    <flowcell
    ID>
    Characters
    allowed:
    a-z, A-Z, 0-9
    4    <lane> Numerical Lane number
    5    <tile> Numerical Tile number
    6    <x_pos> Numerical X coordinate of cluster
    7    <y_pos> Numerical Y coordinate of cluster
    */
    //    A00953:544:HMTFHDSX3:2:1101:6768:1
    //             1        :  2   :    3       :   4  :  5   :   6   :  7
    //    "^([a-zA-Z0-9_-]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:[0-9]+:[0-9-]+:[0-9-]+)$"
    static final Pattern CASAVA_LIKE_EXCLUDE_REGEXP = Pattern.compile(
            "^([a-zA-Z0-9_-]+:[0-9]+:[a-zA-Z0-9_-]+:[0-9]+:[0-9]+:[0-9-]+:[0-9-]+)$");
    // Provided readname structure is @{readkey}{separator:1(.|/|:|_)}{index:1(0:1:2)}
    static final Pattern SPLIT_REGEXP = Pattern.compile("^(.*)(?:[\\.|:|/|_])([0-9]+)$");

    private final String name;
    private final String nameWithoutIndex;
    private final String index;
    private final String bases;
    private final String qualityScores;

    public FastqRead(String name, String bases, String qualityScores) {
        this.name = name;
        this.bases = bases;
        this.qualityScores = qualityScores;

        nameWithoutIndex = getReadPart(name, ReadNameGroup.NAME);
        index = getReadPart(name, ReadNameGroup.INDEX);
    }

    public String getName() {
        return name;
    }

    public String getBases() {
        return bases;
    }

    public String getQualityScores() {
        return qualityScores;
    }

    public String getNameWithoutIndex() {
        return nameWithoutIndex;
    }

    public String getIndex() {
        return index;
    }

    private String
    getReadPart(String readName, ReadNameGroup group) throws ReadWriterException {
        Matcher casavaLikeMatcher = CASAVA_LIKE_EXCLUDE_REGEXP.matcher(readName);
        if (!casavaLikeMatcher.find()) {
            Matcher m = SPLIT_REGEXP.matcher(readName);
            if (m.find()) {
                return m.group(group.getValue());
            }
        }
        if (group == ReadNameGroup.NAME) {
            return readName;
        } else {
            return null;
        }
    }

    public enum ReadNameGroup {
        NAME(1),
        INDEX(2);

        private final int value;

        ReadNameGroup(int value) {
            this.value = value;
        }

        public int getValue() {
            return this.value;
        }
    };
}
