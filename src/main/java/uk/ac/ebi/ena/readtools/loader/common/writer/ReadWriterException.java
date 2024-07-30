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
package uk.ac.ebi.ena.readtools.loader.common.writer;

public class ReadWriterException extends RuntimeException {
  private static final long serialVersionUID = 1L;
  private final ErrorType errorType;

  public ReadWriterException(String value, ErrorType errorType) {
    super(value);
    this.errorType = errorType;
  }

  public ReadWriterException(ErrorType errorType) {
    super();
    this.errorType = errorType;
  }

  public ReadWriterException(Throwable cause, ErrorType errorType) {
    super(cause);
    this.errorType = errorType;
  }

  public ErrorType getErrorType() {
    return errorType;
  }

  public enum ErrorType {
    UNEXPECTED_PAIR_NUMBER,
    SPOT_DUPLICATE,
    SORTING_ERROR,
    BASES_QUALITIES_LENGTH_MISMATCH,
    SAM_RECORD_ERROR,
    INVALID_READ_NAME
  }
}
