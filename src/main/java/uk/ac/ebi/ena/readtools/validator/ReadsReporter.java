package uk.ac.ebi.ena.readtools.validator;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.List;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.ScannerMessage;
import uk.ac.ebi.ena.readtools.webin.cli.rawreads.ScannerMessage.ScannerErrorMessage;

public class ReadsReporter
{
  public void write(File file, Severity severity, String origin, String message) {
    try (PrintWriter writer = new PrintWriter(file)) {
      String out = String.format("%s: %s%s", severity.toString(), origin, message);
      writer.write(out);
      writer.flush();
    } catch (FileNotFoundException e) {
      System.out.println(e);
    }
  }

  public void write(File file, ScannerMessage scannerMessage) {
    Severity severity = scannerMessage instanceof ScannerErrorMessage ? Severity.ERROR : Severity.INFO;
    String msg = scannerMessage.getMessage();
    if (null != scannerMessage.getThrowable())
      msg += scannerMessage.getThrowable().toString();
    write(file, severity, scannerMessage.getOrigin(), msg);
  }

  public void write(File file, List<ScannerMessage> messages) {
    for (ScannerMessage sm: messages)
      write(file, sm);
  }

  public enum Severity {
    ERROR(3),
    WARNING(2),
    INFO(1),
    FIX(0);

    private Integer intVal;

    Severity(Integer intVal) {
      this.intVal = intVal;
    }

    public Integer getIntVal() {
      return intVal;
    }
  }
}
