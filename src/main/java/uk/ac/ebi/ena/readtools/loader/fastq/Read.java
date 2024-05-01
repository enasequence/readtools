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
package uk.ac.ebi.ena.readtools.loader.fastq;

import java.io.Serializable;
import uk.ac.ebi.ena.readtools.loader.common.writer.Spot;

/** Holds raw unpaired read information. */
public class Read implements Serializable, Spot {
  protected final String name;
  protected final String bases;
  protected final String qualityScores;
  protected final String defaultReadIndex;

  public Read(String name, String bases, String qualityScores) {
    this(name, bases, qualityScores, null);
  }

  public Read(String name, String bases, String qualityScores, String defaultReadIndex) {
    this.name = name;
    this.bases = bases;
    this.qualityScores = qualityScores;
    this.defaultReadIndex = defaultReadIndex;
  }

  @Override
  public String getName() {
    return name;
  }

  @Override
  public long getBaseCount() {
    return bases.length();
  }

  public String getBases() {
    return bases;
  }

  public String getQualityScores() {
    return qualityScores;
  }

  public String getDefaultReadIndex() {
    return defaultReadIndex;
  }

  public String toString() {
    return "base_name = ["
        + name
        + "]\n"
        + "bases = ["
        + bases
        + "], length = "
        + (null == bases ? "null" : bases.length())
        + "\nquals = ["
        + qualityScores
        + "], length = "
        + (null == qualityScores ? "null" : qualityScores.length());
  }

  public long getSizeBytes() {
    return name.getBytes().length + bases.getBytes().length + qualityScores.getBytes().length;
  }
}
