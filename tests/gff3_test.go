package main

import (
	"os"
	"strings"
	"testing"

	"../../gff3"
)

func TestSimpleParseLineWithReader(t *testing.T) {
	fakeGff3Line := "chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tID=ENSG00000223972.5;gene_id=ENSG00000223972.5;gene_type=transcribed_unprocessed_pseudogene;gene_status=KNOWN;gene_name=DDX11L1;level=2;havana_gene=OTTHUMG00000000961.2"
	myReader := gff3.NewReader(strings.NewReader(fakeGff3Line))
	if myRecord, err := myReader.Read(); err != nil {
		t.Errorf("record was not correctly parsed and returned an error")
	} else if myRecord.Complete != true {
		t.Errorf("record was not correctly parsed, but did not throw an error")
	}
}

func TestMain(m *testing.M) {
	os.Exit(m.Run())
}
