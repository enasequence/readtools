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
package uk.ac.ebi.ena.readtools.sam;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.zip.GZIPOutputStream;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.FileConverter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.ref.CRAMReferenceSource;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;

import uk.ac.ebi.ena.readtools.cram.CramTools;
import uk.ac.ebi.ena.readtools.cram.ref.ENAReferenceSource;
import uk.ac.ebi.ena.readtools.fastq.MultiFastqOutputter;

public class Sam2Fastq {
	private static final Log log = Log.getInstance(Sam2Fastq.class);

	public static final boolean INCLUDE_NON_PRIMARY_ALIGNMENTS = false;
	public static final boolean INCLUDE_NON_PF_READS = false;

	private static void printUsage(JCommander jc) {
		StringBuilder sb = new StringBuilder();
		sb.append("\n");
		jc.usage(sb);

		System.out.println("Version " + Sam2Fastq.class.getPackage().getImplementationVersion());
		System.out.println(sb.toString());
	}
	
	public static void main(String[] args) throws Exception {
		Params params = new Params();
		JCommander jc = new JCommander(params);
		try {
			jc.parse(args);
		} catch (Exception e) {
			System.out.println("Failed to parse parameteres, detailed message below: ");
			System.out.println(e.getMessage());
			System.out.println();
			System.out.println("See usage: -h");
			System.exit(1);
		}

		if (args.length == 0 || params.help) {
			printUsage(jc);
			System.exit(1);
		}

		Log.setGlobalLogLevel( params.logLevel );

		new Sam2Fastq().create(params);
	}

	public void create( Params params ) throws Exception {
		CRAMReferenceSource referenceSource = new ENAReferenceSource( /* params.reference */ );

		if( params.reference == null )
			log.warn( "No reference file specified, remote access over internet may be used to download public sequences. " );

		final AtomicBoolean brokenPipe = new AtomicBoolean( false );
		try
		{
			sun.misc.Signal.handle( new sun.misc.Signal( "PIPE" ), new sun.misc.SignalHandler()
			{
				@Override public void
				handle( sun.misc.Signal sig )
				{
					brokenPipe.set( true );
				}
			} );
		}catch ( IllegalArgumentException e )
		{
			e.printStackTrace();
		}

		CollatingDumper d = new CollatingDumper( params.samFile, referenceSource, params.nofStreams, params.fastqBaseName, params.gzip,
				params.maxRecords, params.reverse, params.defaultQS, brokenPipe );
		d.prefix = params.prefix;
		d.run();

		if( d.exception != null )
			throw d.exception;
	}
	
	private static abstract class Dumper implements Runnable {
		protected File samFile;
		protected byte[] ref = null;
		protected CRAMReferenceSource referenceSource;
		protected FileOutput[] outputs;
		protected long maxRecords = -1;
		protected SAMFileHeader samHeader;
		protected Exception exception;
		private boolean reverse = false;
		protected AtomicBoolean brokenPipe;

		public Dumper(File samFile, CRAMReferenceSource referenceSource, int nofStreams, String fastqBaseName,
				boolean gzip, long maxRecords, boolean reverse, int defaultQS, AtomicBoolean brokenPipe)
				throws IOException {

			this.samFile = samFile;
			this.referenceSource = referenceSource;
			this.maxRecords = maxRecords;
			this.reverse = reverse;
			this.brokenPipe = brokenPipe;
			outputs = new FileOutput[nofStreams];
			for (int index = 0; index < outputs.length; index++)
				outputs[index] = new FileOutput();

			if (fastqBaseName == null) {
				OutputStream joinedOS = System.out;
				if (gzip)
					joinedOS = (new GZIPOutputStream(joinedOS));
				for (int index = 0; index < outputs.length; index++)
					outputs[index].outputStream = joinedOS;
			} else {
				String extension = ".fastq" + (gzip ? ".gz" : "");
				String path;
				for (int index = 0; index < outputs.length; index++) {
					if (index == 0)
						path = fastqBaseName + extension;
					else
						path = fastqBaseName + "_" + index + extension;

					outputs[index].file = new File(path);
					OutputStream os = new BufferedOutputStream(new FileOutputStream(outputs[index].file));

					if (gzip)
						os = new GZIPOutputStream(os);

					outputs[index].outputStream = os;
				}
			}
		}

		protected abstract MultiFastqOutputter createFastqWriter();

		protected void doRun() throws IOException {
			SamReaderFactory.setDefaultValidationStringency(ValidationStringency.LENIENT);

			final SamReader samReader = SamReaderFactory.makeDefault()
					.referenceSource(referenceSource)
					.open(samFile);

			samHeader = samReader.getFileHeader();

			MultiFastqOutputter fastqOutputter = createFastqWriter();

			for (final SAMRecord currentRecord : samReader) {

				if (currentRecord.isSecondaryOrSupplementary() && !INCLUDE_NON_PRIMARY_ALIGNMENTS) {
					continue;
				}

				// Skip non-PF reads as necessary
				if (currentRecord.getReadFailsVendorQualityCheckFlag() && !INCLUDE_NON_PF_READS) {
					continue;
				}

				String readName = currentRecord.getReadName();

				byte[] readBases = Arrays.copyOf(currentRecord.getReadBases(), currentRecord.getReadBases().length);
				byte[] baseQualities = currentRecord.getBaseQualityString().getBytes(StandardCharsets.UTF_8);

				if (reverse && currentRecord.getReadNegativeStrandFlag()) {
					SequenceUtil.reverseComplement(readBases);
					SequenceUtil.reverseQualities(baseQualities);
				}

				fastqOutputter.writeRead(
						readName.getBytes(StandardCharsets.UTF_8), currentRecord.getFlags(), readBases, baseQualities);
			}

			CloserUtil.close(samReader);

			if (!brokenPipe.get())
				fastqOutputter.finish();
		}

		@Override
		public void run() {
			try {
				doRun();

				if (outputs != null) {
					for (FileOutput os : outputs)
						os.close();
				}
			} catch (Exception e) {
				this.exception = e;
			}
		}
	}
	
	private static class CollatingDumper extends Dumper {
		private FileOutput fo = new FileOutput();
		private String prefix;
		private long counter = 1;
		private MultiFastqOutputter multiFastqOutputter;
		private int defaultQS;

		public CollatingDumper(File samFile, CRAMReferenceSource referenceSource, int nofStreams,
				String fastqBaseName, boolean gzip, long maxRecords, boolean reverse, int defaultQS,
				AtomicBoolean brokenPipe) throws IOException {
			super(samFile, referenceSource, nofStreams, fastqBaseName, gzip, maxRecords, reverse, defaultQS, brokenPipe);
			this.defaultQS = defaultQS;
			this.brokenPipe = brokenPipe;
			fo.file = File.createTempFile(fastqBaseName == null ? "overflow.bam" : fastqBaseName + ".overflow.bam",
					".tmp");
			fo.file.deleteOnExit();
			fo.outputStream = new BufferedOutputStream(new FileOutputStream(fo.file));
		}

		@Override
		protected MultiFastqOutputter createFastqWriter() {
			multiFastqOutputter = new MultiFastqOutputter(outputs, fo, samHeader);
			if (prefix != null) {
				multiFastqOutputter.setPrefix(prefix.getBytes());
			}
			return multiFastqOutputter;
		}

		@Override
		public void doRun() throws IOException {
			super.doRun();

			fo.close();

			if (fo.empty)
				return;

			log.info("Sorting overflow BAM: ", fo.file.length());
			SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
			SamReader r = SamReaderFactory.makeDefault().open(fo.file);
			SAMRecordIterator iterator = r.iterator();
			if (!iterator.hasNext()) {
				r.close();
				fo.file.delete();
				return;
			}

			SAMRecord r1 = iterator.next();
			SAMRecord r2 = null;
			counter = multiFastqOutputter.getCounter();
			log.info("Counter=" + counter);
			while (!brokenPipe.get() && iterator.hasNext()) {
				r2 = iterator.next();
				if (r1.getReadName().equals(r2.getReadName())) {
					print(r1, r2);
					counter++;
					r1 = null;
					if (!iterator.hasNext())
						break;
					r1 = iterator.next();
					r2 = null;
				} else {
					print(r1, 0);
					r1 = r2;
					r2 = null;
					counter++;
				}
			}
			if (r1 != null)
				print(r1, 0);
			r.close();
			fo.file.delete();
		}

		private void print(SAMRecord r1, SAMRecord r2) throws IOException {
			if (r1.getFirstOfPairFlag()) {
				print(r1, 1);
				print(r2, 2);
			} else {
				print(r1, 2);
				print(r2, 1);
			}
		}

		private void print(SAMRecord r, int index) throws IOException {
			OutputStream os = outputs[index];
			os.write('@');
			if (prefix != null) {
				os.write(prefix.getBytes());
				os.write('.');
				os.write(String.valueOf(counter).getBytes());
				os.write(' ');
			}
			os.write(r.getReadName().getBytes());
			if (index > 0) {
				os.write('/');
				os.write(48 + index);
			}
			os.write('\n');
			os.write(r.getReadBases());
			os.write("\n+\n".getBytes());
			os.write(r.getBaseQualityString().getBytes());
			os.write('\n');
		}
	}

	private static class FileOutput extends OutputStream {
		File file;
		OutputStream outputStream;
		boolean empty = true;

		@Override
		public void write(int b) throws IOException {
			outputStream.write(b);
			empty = false;
		}

		@Override
		public void write(byte[] b) throws IOException {
			outputStream.write(b);
			empty = false;
		}

		@Override
		public void write(byte[] b, int off, int len) throws IOException {
			outputStream.write(b, off, len);
			empty = false;
		}

		@Override
		public void flush() throws IOException {
			outputStream.flush();
		}

		@Override
		public void close() throws IOException {
			if (outputStream != null && outputStream != System.out && outputStream != System.err) {
				outputStream.close();
				outputStream = null;
			}
			if (empty && file != null && file.exists())
				file.delete();
		}
	}

	@Parameters(commandDescription = "SAM to FastQ dump conversion. ")
	public static class Params {
		@Parameter(names = { "-l", "--log-level" }, description = "Change log level: DEBUG, INFO, WARNING, ERROR.", converter = CramTools.LevelConverter.class)
		Log.LogLevel logLevel = Log.LogLevel.ERROR;

		@Parameter(names = { "-h", "--help" }, description = "Print help and quit")
		boolean help = false;

		@Parameter(names = { "--input-sam-file", "-I" }, converter = FileConverter.class, description = "The path to the SAM file to uncompress. Omit if standard input (pipe).")
		File samFile;

//		@Parameter(names = { "--reference-fasta-file", "-R" }, converter = FileConverter.class, description = "Path to the reference fasta file, it must be uncompressed and indexed (use 'samtools faidx' for example). ")
		File reference;

		@Parameter(names = { "--fastq-base-name", "-F" }, description = "'_number.fastq[.gz] will be appended to this string to obtain output fastq file name. If this parameter is omitted then all reads are printed with no garanteed order.")
		String fastqBaseName;

		@Parameter(names = { "--gzip", "-z" }, description = "Compress fastq files with gzip.")
		boolean gzip;

		@Parameter(names = { "--reverse" }, description = "Re-reverse reads mapped to negative strand.")
		boolean reverse = false;

		@Parameter(names = { "--enumerate" }, description = "Append read names with read index (/1 for first in pair, /2 for second in pair).")
		boolean appendSegmentIndexToReadNames;

		@Parameter(names = { "--max-records" }, description = "Stop after reading this many records.")
		long maxRecords = -1;

		@Parameter(names = { "--read-name-prefix" }, description = "Replace read names with this prefix and a sequential integer.")
		String prefix = null;

		@Parameter(names = { "--default-quality-score" }, description = "Use this quality score (decimal representation of ASCII symbol) as a default value when the original quality score was lost due to compression. Minimum is 33.")
		int defaultQS = '?';

		@Parameter(names = { "--ignore-md5-mismatch" }, description = "Issue a warning on sequence MD5 mismatch and continue. This does not garantee the data will be read succesfully. ")
		public boolean ignoreMD5Mismatch = false;

		@Parameter(names = { "--skip-md5-check" }, description = "Skip MD5 checks when reading the header.")
		public boolean skipMD5Checks = false;
		
		public int nofStreams = 3;
	}
}
