package main

import (
	"os"
	"strings"
	"testing"

	"../../gff3"
)

var fakeGff3Line = "chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tID=ENSG00000223972.5;gene_id=ENSG00000223972.5;gene_type=transcribed_unprocessed_pseudogene;gene_status=KNOWN;gene_name=DDX11L1;level=2;havana_gene=OTTHUMG00000000961.2"

func TestSimpleParseLineWithReader(t *testing.T) {
	myReader := gff3.NewReader(strings.NewReader(fakeGff3Line))
	if myRecord, err := myReader.Read(); err != nil {
		t.Errorf("record was not correctly parsed and returned an error")
	} else if !myRecord.Complete {
		t.Errorf("record was not correctly parsed, but did not throw an error")
	}
}

func TestParseTwoLinesWithReader(t *testing.T) {
	fakeGff3Line := "chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tID=ENSG00000223972.5;gene_id=ENSG00000223972.5;gene_type=transcribed_unprocessed_pseudogene;gene_status=KNOWN;gene_name=DDX11L1;level=2;havana_gene=OTTHUMG00000000961.2\nchr1\tHAVANA\tgene\t69091\t70008\t.\t+\t.\tID=ENSG00000186092.4;gene_id=ENSG00000186092.4;gene_type=protein_coding;gene_status=KNOWN;gene_name=OR4F5;level=2;havana_gene=OTTHUMG00000001094.2"
	myReader := gff3.NewReader(strings.NewReader(fakeGff3Line))
	if myRecord, err := myReader.Read(); err != nil {
		t.Errorf("first record was not correctly parsed and returned an error")
	} else if !myRecord.Complete {
		t.Errorf("first record was not correctly parsed, but did not throw an error")
	}
	if myRecord, err := myReader.Read(); err != nil {
		t.Errorf("second record was not correctly parsed and returned an error")
	} else if !myRecord.Complete {
		t.Errorf("second record was not correctly parsed, but did not throw an error")
	}
}

func TestParseSeveralCommentsBeforeLine(t *testing.T) {
	fakeGff3Line := "#GFF3 file\n#header information\n#file version\nchr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tID=ENSG00000223972.5;gene_id=ENSG00000223972.5;gene_type=transcribed_unprocessed_pseudogene;gene_status=KNOWN;gene_name=DDX11L1;level=2;havana_gene=OTTHUMG00000000961.2"
	myReader := gff3.NewReader(strings.NewReader(fakeGff3Line))
	if myRecord, err := myReader.Read(); err != nil {
		t.Errorf("record was not correctly parsed and returned an error")
	} else if !myRecord.Complete {
		t.Errorf("record was not correctly parsed, but did not throw an error")
	}
}

func TestRecordFilterStrand(t *testing.T) {
	myReader := gff3.NewReader(strings.NewReader(fakeGff3Line))
	myRecord, _ := myReader.Read()
	filtRecord := myRecord.FilterByField("strand", "+")
	if !filtRecord.Complete {
		t.Errorf("record did not successfully pass filter")
	}
	filtRecord = myRecord.FilterByField("strand", "-")
	if filtRecord.Complete {
		t.Errorf("record did not appropriately fail filter")
	}
}

func TestRecordFilterAttribute(t *testing.T) {
	myReader := gff3.NewReader(strings.NewReader(fakeGff3Line))
	myRecord, _ := myReader.Read()
	filtRecord := myRecord.FilterByAttribute("level", "2")
	if !filtRecord.Complete {
		t.Errorf("record did not successfully pass filter")
	}
	filtRecord = myRecord.FilterByAttribute("gene", "GAPDH")
	if filtRecord.Complete {
		t.Errorf("record did not appropriately fail fitler")
	}
}

func TestChainedRecordFilters(t *testing.T) {
	myReader := gff3.NewReader(strings.NewReader(fakeGff3Line))
	myRecord, _ := myReader.Read()
	filtRecord := myRecord.FilterByAttribute("level", "2").FilterByField("strand", "+").FilterByField("type", "gene")
	if !filtRecord.Complete {
		t.Errorf("record did not successfully pass filter")
	}
	filtRecord = myRecord.FilterByAttribute("level", "2").FilterByField("strand", "+").FilterByField("type", "exon")
	if filtRecord.Complete {
		t.Errorf("record did not appropriately fail filter")
	}
	filtRecord = myRecord.FilterByAttribute("level", "2").FilterByField("strand", "-").FilterByField("type", "gene")
	if filtRecord.Complete {
		t.Errorf("record did not appropriately fail filter")
	}
}

func TestMain(m *testing.M) {
	os.Exit(m.Run())
}
