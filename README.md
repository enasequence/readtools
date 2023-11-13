# Readtools

Readtools is a Java library for reading and validating read data file formats including Fastq, BAM and CRAM. It is used by Webin-CLI and ENA's internal processing pipelines.

# License

Copyright 2015-2023 EMBL - European Bioinformatics Institute Licensed under the Apache License, Version 2.0 (the "License"); 
you may not use this file except in compliance with the License. 
You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0 
Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

| format | parser           | feeder                            | common validation                                                                     |specific validation|
|--------|------------------|-----------------------------------|---------------------------------------------------------------------------------------|---|
| fastq  | ReadReader       | AutoNormalizeQualityReadConverter |                                                                                       | read-qual len, qality format, read name format, pairing, duplicates checks            |
| sam    | htsjdk.SamReader | direct                            | 1+ reads, no empty reads, IUPAC nucleotides, qality level, all or no bases have quals | header, references, read-qual len or read-level qual checks |
| fasta  | none             | none                              | reads qual control filtering                                                          |

