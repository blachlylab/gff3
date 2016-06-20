// http://www.sequenceontology.org/gff3.shtml
package gff3

// Imports
// standard library
import (
	"log"
	"strconv"
	"strings"
)

type Record struct {
	Complete        bool // 'true' once all fields loaded
	seqidField      string
	sourceField     string
	typeField       string
	startField      uint
	endField        uint
	scoreField      float64
	strandField     byte
	phaseField      int
	attributesField map[string]string
}

// loads a GFF3 record from a line of input into the Record struct
func LoadRecord(line string) *Record {
	// ignore comment lines
	if line[0] == '#' {
		return new(Record)
	}

	// strip trailing newline, if any
	line = strings.TrimSuffix(line, "\n")

	// split into fields
	fields := strings.Split(line, "\t")
	if len(fields) != 9 {
		// comment lines should already have been dealt with,
		// so this is a malformed record
		log.Fatalln("Malformed record: ", line)
	}

	r := new(Record)
	r.seqidField = fields[0]
	r.sourceField = fields[1]
	r.typeField = fields[2]
	r.startField, _ = strconv.Atoi(fields[3])
	r.endField, _ = strconv.Atoi(fields[4])
	r.scoreField, _ = strconv.ParseFloat(fields[5], 64)
	r.strandField = fields[6][0] // one byte char: +, -, ., or ?
	r.phaseField, _ = strconv.Atoi(fields[7])
	r.attributesField = fields[8]

	// validate is currently stub function always true
	if r.validate() {
		r.Complete = true
		return r
	} else {
		// need to return an error also
		return r
	}
}

// checks validity of GFF3 record
// stub function always true
func (r *Record) validate() bool {
	// validate but do not rely on "Complete" field
	// to do: can I modify/update complete from here since I am being passed a pointer?
	return true
}

// Filter functions
// Field filters and Attribute filters

// FilterByField wraps individual functions that filter on
// various fields of the Record type
// filters are of the "equals exactly" type. They are case sensitive
// presently there are filters for:
// typeField (e.g. exon, gene, transcript)
// strandField (e.g. "+", "-", ".", "?")
func (r *Record) FilterByField(field, value string) *Record {
	// if record is incomplete, nothing to do
	if !r.Complete {
		return r
	}
	switch field {
	case "type", "Type", "typeField", "TypeField":
		return r.filterByTypeField(value)
	case "strand", "Strand", "strandField", "StrandField":
		return r.filterByStrandField(value[0])
	}
	panic("Unimplemented filter field")
	// should not be reached
	return r
}

func (r *Record) filterByTypeField(filterValue string) *Record {
	// not exported, so r.Complete should have already been checked
	if r.typeField == filterValue {
		return r
	} else {
		r.Complete = false
		return r
	}
}

func (r *Record) filterByStrandField(filterValue byte) *Record {
	// not exported, so r.Complete should have already been checked
	if r.strandField == filterValue {
		return r
	} else {
		r.Complete = false
		return r
	}
}

// Attribute filter
// Attributes are stored in the GFF3 file of the form
// key=value or key=value1,value2,... and are semicolon-delimited, e.g.:
// key1=A;key2=B,C;key9=Z
// In the record struct, they are stored in a map key->value(s)
func (r *Record) FilterByAttribute(attribute, filterValue string) *Record {
	// if record is incomplete, nothing to do
	if !r.Complete {
		return r
	}
	// if attribute not found, return incomplete or empty record
	if rawv, ok := r.attributesField[attribute]; !ok {
		r.Complete = false
		return r
	}
	// attribute was found in the map and its value(s) stored in v
	// gff3 spec allows multiple values for each attribute to be
	// separated by commas, e.g.:
	// tag=appris,basic,CCDS
	valSlice := strings.Split(rawv, ",")
	for _, val := range valSlice {
		if val == filterValue {
			return r
		}
	}
	// filterValue was not found among attribute's values
	// return incomplete or empty record
	r.Complete = false
	return r
}
