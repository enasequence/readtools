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
package uk.ac.ebi.ena.readtools.common.reads;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Shared constants and helpers for Casava 1.8 read name parsing.
 *
 * <p>Casava 1.8 header format:
 *
 * <pre>@instrument:run:flowcell:lane:tile:x:y readnum:filter:control:barcode</pre>
 *
 * <p>Regex groups: (1) base name, (2) whitespace, (3) read number, (4) rest (empty or :barcode...).
 */
public final class CasavaRead {

  private static final String CASAVA_18_CORE =
      "(.+)( +|\\t+)([0-9]+):([YN]):([0-9]*[02468])($|:.*$)";

  /** Pattern for processed read names (no @ prefix). Used by validators. */
  public static final Pattern P_CASAVA_18_NAME = Pattern.compile("^" + CASAVA_18_CORE);

  /** Pattern for raw FASTQ lines (with @ prefix). Used by ReadReader. */
  public static final Pattern P_CASAVA_18_RAW_LINE = Pattern.compile("^@" + CASAVA_18_CORE);

  private CasavaRead() {}

  /** Returns group(1) — the base/instrument name — or null if not Casava 1.8 format. */
  public static String getBaseNameOrNull(String readName) {
    Matcher m = P_CASAVA_18_NAME.matcher(readName);
    if (!m.matches()) return null;
    return m.group(1);
  }

  /** Returns group(3) — the read number (1, 2, ...) — or null if not Casava 1.8 format. */
  public static String getReadIndexOrNull(String readName) {
    Matcher m = P_CASAVA_18_NAME.matcher(readName);
    if (!m.matches()) return null;
    return m.group(3);
  }

  /**
   * Returns the barcode/index sequence from the Casava tail, or null if not Casava 1.8 or no
   * barcode present. For example, from "inst:1:FC:1:1:0:0 1:Y:18:ATCACG+GCGCTA" returns
   * "ATCACG+GCGCTA".
   */
  public static String getBarcodeOrNull(String readName) {
    Matcher m = P_CASAVA_18_NAME.matcher(readName);
    if (!m.matches()) return null;
    String tail = m.group(6); // ($|:.*$)
    if (tail.isEmpty()) return null;
    // tail starts with ':', strip it
    String barcode = tail.substring(1);
    if (barcode.isEmpty()) return null;
    return barcode;
  }

  /** Returns true if the filter flag is Y (read failed filtering). */
  public static boolean isFiltered(String readName) {
    Matcher m = P_CASAVA_18_NAME.matcher(readName);
    if (!m.matches()) return false;
    return "Y".equals(m.group(4));
  }
}
