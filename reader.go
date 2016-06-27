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
	ok := r.s.Scan() // read a line, checking for end of file
	// is the scanner empty?
	if !ok {
		return nil, io.EOF
	}
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
	rec.SeqidField = fields[0]
	rec.SourceField = fields[1]
	rec.TypeField = fields[2]
	rec.StartField, _ = strconv.Atoi(fields[3])
	rec.EndField, _ = strconv.Atoi(fields[4])
	rec.ScoreField, _ = strconv.ParseFloat(fields[5], 64)
	rec.StrandField = fields[6][0] // one byte char: +, -, ., or ?
	rec.PhaseField, _ = strconv.Atoi(fields[7])
	// must initialize the nil map or face a runtime panic
	// maybe better to do this here rather than the Record definition
	//   just in case somebody wants to implement Record without the attribute field; it's slow
	rec.AttributesField = make(map[string]string)
	attrs := fields[8]
	var eqIndex int
	for i := strings.Index(attrs, ";"); i > 0; i = strings.Index(attrs, ";") {
		eqIndex = strings.Index(attrs[:i], "=")
		rec.AttributesField[attrs[:i][:eqIndex]] = attrs[:i][eqIndex+1:]
		attrs = attrs[i+1:]
	}
	// no trailing semicolon on gtf lines; grab the last entry
	eqIndex = strings.Index(attrs, "=")
	rec.AttributesField[attrs[:eqIndex]] = attrs[eqIndex+1:]

	// validate is currently stub function always true
	if rec.Validate() {
		rec.Complete = true
		return rec, nil
	} else {
		// need to return an error also
		return rec, nil
	}
}
