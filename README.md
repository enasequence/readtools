# readtools
## ENA read tools

* Typical cramtools usage in processing:

    cram dump: java -XX:+UseSerialGC -Xmx10G -Dsamjdk.use_cram_ref_download=true -Djava.io.tmpdir="/var/tmp" -jar cramtools.jar fastq --reverse --gzip $@
    cram stats: java -XX:+UseSerialGC -Xmx4G -cp cramtools.jar net.sf.cram.CramTools bam -c -F 2304 -I $@
