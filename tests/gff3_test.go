package unit_tests

import (
	"os"
	"reflect"
	"strings"
	"testing"

	"../../gff3"
)

var fakeGff3Line = "chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tID=ENSG00000223972.5;gene_id=ENSG00000223972.5;gene_type=transcribed_unprocessed_pseudogene;gene_status=KNOWN;gene_name=DDX11L1;level=2;havana_gene=OTTHUMG00000000961.2"

func TestSimpleParseLineWithReader(t *testing.T) {
	myReader := gff3.NewReader(strings.NewReader(fakeGff3Line))
	myRecord, err := myReader.Read()
	if err != nil {
		t.Errorf("record was not correctly parsed and returned an error")
	}
	// file opened fine, but is the content correct?
	var fakeGff3Solution = gff3.Record{
		Complete:    true,
		SeqidField:  "chr1",
		SourceField: "HAVANA",
		TypeField:   "gene",
		StartField:  11869,
		EndField:    14409,
		ScoreField:  0,
		StrandField: '+',
		PhaseField:  0,
		AttributesField: map[string]string{
			"ID":          "ENSG00000223972.5",
			"gene_id":     "ENSG00000223972.5",
			"gene_type":   "transcribed_unprocessed_pseudogene",
			"gene_status": "KNOWN",
			"gene_name":   "DDX11L1",
			"level":       "2",
			"havana_gene": "OTTHUMG00000000961.2",
		},
	}
	// deep equal will recursively compare all attributes and values of the structs
	// a simple equal "==" doesn't suffice because keys in the attributes map may be reordered
	if !reflect.DeepEqual(myRecord, &fakeGff3Solution) {
		t.Errorf("record was not correctly parsed, but did not throw an error")
	}
	fakeGff3Solution.StartField = 1
	if reflect.DeepEqual(myRecord, &fakeGff3Solution) {
		t.Errorf("records are different and should not be equal")
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
