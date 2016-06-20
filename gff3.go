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

// A Reader reads records from a GFF3-formatted file
//
// As returned by NewReader, a Reader expects input to conform to
// http://www.sequenceontology.org/gff3.shtml
//
// As presently (2016-06-20) implemented, there may be some
// GENCODE or Ensembl specific attributes expected hereunder
//
// TODO: copy similar boilerplate from golang.org/src/encoding/csv/reader.go
//
type Reader struct {
	s	*bufio.Scanner
	line	int
	column	int
	field	bytes.Buffer
}

// NewReader returns a new Reader that reads from r
// This is simply a wrapper to buffer input io.Reader if not already buffered
// also use the Scanner which is apparently idiomatic way to read lines now
func NewReader(r io.Reader) *Reader {
	return &Reader{
		s: bufio.NewScanner(r), // buffer input
	}
}

func (r *Reader) Read() (record *Record, err error) {
	// parseRecord one at a time, but retry if nil record (e.g. comment)
	for {
		record, err = r.parseRecord()
		// nil record, nil error could represent a comment line
		if record != nil {
			break
		}
		if err != nil {
			return nil, err
		}
	}
	return record, nil
}

// func (r *Reader) ReadAll() (records []Record, err error) 

// parseRecord parses the elements of a single row into a GFF3 struct
func (r *Reader) parseRecord() (*Record, err error) {
	// line numbering starts with 1, not 0, so increment straightaway
	r.line++
	_ = r.s.Scan()	// read a line

	// Examine first byte -- if comment, skip this row
	if r.s.Bytes()[0] == '#' {
		return nil, nil	// no record, no error
		// returning nil might break Filtering which expects *Record
		// however I think the for {} loop in Read() will keep looking
		// for a valid record
	}

	// I don't think I need to strip trailing newline

	// split into fields
	// this will allocate a new string and a new array (2 alloc/op)
	fields := strings.Split(r.s.Text(), "\t")
	if len(fields != 9) {
		// comment lines should have already been dealt with,
		// so this is a malformed record
		log.Fatalln("Malformed record: ", r.s.Text())
	}

	rec := new(Record)
	rec.seqidField = fields[0]
	rec.sourceField = fields[1]
	rec.typeField = fields[2]
	rec.startField, _ = strconv.Atoi(fields[3])
	rec.endField, _ = strconv.Atoi(fields[4])
	rec.scoreField, _ = strconv.ParseFloat(fields[5], 64)
	rec.strandField = fields[6][0] // one byte char: +, -, ., or ?
	rec.phaseField, _ = strconv.Atoi(fields[7])
	rec.attributesField = fields[8]

	// validate is currently stub function always true
	if rec.validate() {
		rec.Complete = true
		return rec, nil
	} else {
		// need to return an error also
		return rec, nil
	}
}

// loads a GFF3 record from a line of input into the Record struct
// Wrote this before conversion to Reader object
// deprecated - prefixed with x to avoid export
func xLoadRecord(line string) *Record {
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
