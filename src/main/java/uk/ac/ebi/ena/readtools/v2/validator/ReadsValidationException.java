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
package uk.ac.ebi.ena.readtools.v2.validator;

public class ReadsValidationException extends Exception {
  private final String errorMessage;
  private final long readIndex;
  private final String readName;

  public ReadsValidationException(String errorMessage) {
    this(errorMessage, 0, "");
  }

  public ReadsValidationException(String errorMessage, long readIndex) {
    this(errorMessage, readIndex, "");
  }

  public ReadsValidationException(String errorMessage, long readIndex, String readName) {
    super(errorMessage);
    this.errorMessage = errorMessage;
    this.readIndex = readIndex;
    this.readName = readName;
  }

  public String getErrorMessage() {
    String message = errorMessage;
    if (readIndex > 0) {
      message = "Read Index: " + readIndex + " " + message;
    }
    if (readName != null) {
      message = "Read Name: " + readName + " " + message;
    }
    return message;
  }

  public long getReadIndex() {
    return readIndex;
  }
}
