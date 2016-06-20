// Package gff3 implements data structure that represents GFF3 formatted
// (http://www.sequenceontology.org/gff3.shtml) data, as well as Reader
// and Writer for this data
//
// Some code inspired by the golang standard library's csv package
package gff3

// Imports
// standard library
import (
	"bufio"
	"bytes"
	"io"
	"log"
	"strconv"
	"strings"
)

// A Record represents a single element in a GFF3 file -- a single row
//
// it can be received from Read or passed to Write
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
