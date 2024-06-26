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
package uk.ac.ebi.ena.readtools.webin.cli.rawreads.refs;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamFiles;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.Log.LogLevel;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.net.URL;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;
import javax.net.ssl.HttpsURLConnection;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class CramReferenceInfo {
  private static final String service_link =
      "https://www.ebi.ac.uk/ena/cram/sequence/%32s/metadata";
  private static final Logger log = LoggerFactory.getLogger(CramReferenceInfo.class);

  public static void main(String[] args) {
    CramReferenceInfo m = new CramReferenceInfo();
    m.fetchReferenceMetadata(args[0] /* "a0d9851da00400dec1098a9255ac712e"*/, System.out);
  }

  public ReferenceInfo fetchReferenceMetadata(String md5, PrintStream... pss) {
    for (int i = 0; i < 4; ++i)
      try {
        String data = fetchData(new URL(String.format(service_link, md5)));
        ReferenceInfo info = parseNCBIJson(data);

        if (null != pss) Arrays.stream(pss).forEachOrdered(ps -> ps.print(data));

        return info;

      } catch (Throwable t) {
        try {
          Thread.sleep(1000);
        } catch (InterruptedException ie) {
          Thread.currentThread().interrupt();
          break;
        }
      }

    return new ReferenceInfo();
  }

  public Map<String, Boolean> confirmFileReferences(File file) throws IOException {
    ThreadPoolExecutor es =
        new ThreadPoolExecutor(10, 10, 1, TimeUnit.HOURS, new ArrayBlockingQueue<>(10));
    Map<String, Boolean> result = new ConcurrentHashMap<>();
    Log.setGlobalLogLevel(LogLevel.ERROR);
    SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
    SamReaderFactory factory = SamReaderFactory.make();
    factory.validationStringency(ValidationStringency.SILENT);
    factory.referenceSource(new ReferenceSource((File) null));
    SamInputResource ir = SamInputResource.of(file);
    File indexMaybe = SamFiles.findIndex(file);
    // System.out.println( "proposed index: " + indexMaybe );
    AtomicLong count = new AtomicLong();
    try (SamReader reader = factory.open(ir)) {
      String msg =
          String.format(
              "Checking reference existence in the CRAM reference registry for %s\n",
              file.getPath());
      log.info(msg);

      es.prestartAllCoreThreads();
      for (SAMSequenceRecord sequenceRecord :
          reader.getFileHeader().getSequenceDictionary().getSequences()) {
        if ("*".equals(sequenceRecord.getSequenceName())) continue;

        es.getQueue()
            .put(
                () -> {
                  count.incrementAndGet();

                  String md5 = sequenceRecord.getAttribute(SAMSequenceRecord.MD5_TAG);
                  if (md5 == null) {
                    throw new IllegalArgumentException(
                        "File header missing MD5 attribute : " + file.getAbsolutePath());
                  }

                  result.computeIfAbsent(
                      md5,
                      k -> {
                        ReferenceInfo info = fetchReferenceMetadata(md5);
                        return k.equals(info.getMd5());
                      });

                  if (0 == count.get() % 10) printProcessedReferenceNumber(count);
                });
      }

      es.shutdown();
      es.awaitTermination(1, TimeUnit.HOURS);

      printProcessedReferenceNumber(count);
      printFlush("\n");

    } catch (InterruptedException e) {
      Thread.currentThread().interrupt();
    }
    return result;
  }

  public void printFlush(String msg) {
    System.out.print(msg);
    System.out.flush();
  }

  public static class ReferenceInfo {
    private String id, md5;

    public String getId() {
      return id;
    }

    public String getMd5() {
      return md5;
    }
  }

  private String fetchData(URL url) throws IOException {
    HttpsURLConnection con = (HttpsURLConnection) url.openConnection();
    StringBuilder sb = new StringBuilder();
    try (BufferedReader br = new BufferedReader(new InputStreamReader(con.getInputStream()))) {
      String input;
      while ((input = br.readLine()) != null) sb.append(input).append('\n');
    }
    return sb.toString();
  }

  private ReferenceInfo parseNCBIJson(String json) {
    ObjectMapper mapper = new ObjectMapper();

    try {
      JsonNode rootNode = mapper.readTree(json);
      JsonNode metadataNode = rootNode.path("metadata");

      Map<String, String> metadata = new HashMap<>();
      metadataNode
          .fieldNames()
          .forEachRemaining((String key) -> metadata.put(key, metadataNode.get(key).asText()));

      ReferenceInfo info = new ReferenceInfo();
      info.id = metadata.get("id");
      info.md5 = metadata.get("md5");

      return info;
    } catch (IOException e) {
      e.printStackTrace();
      throw new RuntimeException(e);
    }
  }

  private void printProcessedReferenceNumber(AtomicLong count) {
    String msg = String.format("\rChecked %16d references(s)", count.get());
    printFlush(msg);
  }
}
