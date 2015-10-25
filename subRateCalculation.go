package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"os"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/linear"
)

func subs(consensus, s seq.Sequence) (int, int) {
	var hamming int
	length := 0

	alpha := s.Alphabet()
	indexOf := alpha.LetterIndex()
	gap := alpha.Gap()
	for c := 0; c < s.Len(); c++ {
		l := s.At(c).L
		if l == gap {
			continue
		}
		if indexOf[l] != indexOf[consensus.At(c).L] {
			hamming++
		}
		length++
	}
	return hamming, length
}

var (
	cLoc string
	sLoc string
)

func main() {
	flag.StringVar(&cLoc, "cLoc", "", "File location of consensus.")
	flag.StringVar(&sLoc, "sLoc", "", "File location of consensus.")
	flag.Parse()

	// Read in sequences
	sFile, err := os.Open(sLoc)
	if err != nil {
		log.Fatal(err)
	}
	defer sFile.Close()

	var ListOfSequences []seq.Sequence
	reader := fasta.NewReader(sFile, linear.NewSeq("", nil, alphabet.DNA))
	sc := seqio.NewScanner(reader)
	for sc.Next() {
		ListOfSequences = append(ListOfSequences, sc.Seq())
	}
	if sc.Error() != nil {
		log.Fatalf("failed to read sequences %v", sc.Error())
	}

	//Read in consensus
	cFile, err := os.Open(cLoc)
	if err != nil {
		log.Fatal(err)
	}
	defer cFile.Close()
	var consensus seq.Sequence
	cReader := fasta.NewReader(cFile, linear.NewSeq("", nil, alphabet.DNA))
	cSc := seqio.NewScanner(cReader)
	for cSc.Next() {
		consensus = cSc.Seq()
		break
	}
	if cSc.Error() != nil {
		log.Fatalf("failed to read consensus: %v", err)
	}

	// Creating a file for the output
	file := fmt.Sprintf("%s_subrate.txt", sLoc)
	out, err := os.Create(file)
	if err != nil {
		log.Fatalf("failed to create %s: %v", file, err)
	}
	defer out.Close()
	fmt.Fprintf(out, "subs\td\n\n")

	for i, b := range ListOfSequences {
		subs, length := subs(consensus, b)
		//calculate JC
		p := float64(subs) / float64(length)
		d := -0.75 * math.Log(1-1.25*p)

		fmt.Fprintf(out, "%v\t%v\n", subs, d)

		fmt.Printf("JC distance %v: %v\n", i+1, d)

	}

}
