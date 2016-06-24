package main

import (
	"os"
    "log"
    
    "../../gff3"
    "github.com/davecheney/profile"
)

func TestRealLifeGFF3Parse(fnIn string) {
	fn, err := os.Open(fnIn)
    if err != nil {
        log.Fatalln("failed to parse", fnIn)
    }
    myReader := gff3.NewReader(fn)
    var myRecord *gff3.Record 
    for {
        myRecord, err = myReader.Read()
        if err != nil {
            break
        }
        if myRecord.Complete {
            continue
        }
    }

}

func main() {
    defer profile.Start(profile.CPUProfile).Stop()
    if len(os.Args) < 2 {
        log.Fatal("failed to supply a gtf. try 'go test sitautional_test.go -args ${gtf}'")
    }
	TestRealLifeGFF3Parse(os.Args[1])
}
