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
	s      *bufio.Scanner
	line   int
	column int
	field  bytes.Buffer
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
func (r *Reader) parseRecord() (rec *Record, err error) {
	// line numbering starts with 1, not 0, so increment straightaway
	r.line++
	_ = r.s.Scan() // read a line

	// Examine first byte -- if comment, skip this row
	if r.s.Bytes()[0] == '#' {
		return nil, nil // no record, no error
		// returning nil might break Filtering which expects *Record
		// however I think the for {} loop in Read() will keep looking
		// for a valid record
	}

	// I don't think I need to strip trailing newline

	// split into fields
	// this will allocate a new string and a new array (2 alloc/op)
	fields := strings.Split(r.s.Text(), "\t")
	if len(fields) != 9 {
		// comment lines should have already been dealt with,
		// so this is a malformed record
		log.Fatalln("Malformed record: ", r.s.Text())
	}

	rec = new(Record)
	rec.seqidField = fields[0]
	rec.sourceField = fields[1]
	rec.typeField = fields[2]
	rec.startField, _ = strconv.Atoi(fields[3])
	rec.endField, _ = strconv.Atoi(fields[4])
	rec.scoreField, _ = strconv.ParseFloat(fields[5], 64)
	rec.strandField = fields[6][0] // one byte char: +, -, ., or ?
	rec.phaseField, _ = strconv.Atoi(fields[7])
	// must initialize the nil map or face a runtime panic
	rec.attributesField = make(map[string]string)
	for _, attribute := range strings.Split(fields[8], ";") {
		kv := strings.Split(attribute, "=")
		key := k[0]
		value := v[1]
		rec.attributesField[key] = value
	}

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
/*
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
	//r.attributesField = fields[8]

	// validate is currently stub function always true
	if r.validate() {
		r.Complete = true
		return r
	} else {
		// need to return an error also
		return r
	}
}
*/
