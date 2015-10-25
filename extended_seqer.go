// Copyright ©2014 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"io"
	"log"
	"math/rand"
	"os"
	"path/filepath"
	//	"strings"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/biogo/seq/multi"
	"github.com/biogo/external/muscle"
)

var (
	cDir      string
	aDir      string
	sDir      string
	maxFam    int
	subSample bool
	consFasta bool
)

func main() {
	flag.IntVar(&maxFam, "maxFam", 0, "maxFam indicates maximum family size considered (0 == no limit).")
	flag.BoolVar(&subSample, "subsample", false, "choose maxFam members of a family if the family has more than maxFam members.")
	flag.BoolVar(&consFasta, "fasta", false, "output consensus as fasta with quality case filtering.")
	flag.StringVar(&cDir, "cDir", "", "target directory for consensus output. If not empty Dir is deleted first.")
	flag.StringVar(&aDir, "aDir", "", "target directory for alignment output. If not empty dir is deleted first.")
	flag.StringVar(&sDir, "sDir", "", "target directory for sequence information output. If not empty dir is deleted first.")
	flag.Parse()

	fmt.Printf("Initialising Files: %s\n", flag.Args()[0])

	checks(cDir, aDir, sDir)

	//Opening files
	f, fErr := os.Open(flag.Args()[0])
	if fErr != nil {
		log.Printf("error: could not open %s to read %v", flag.Args()[0], fErr)
	}
	defer f.Close()

	// Reading in sequences
	var v []seq.Sequence
	r := fasta.NewReader(f, linear.NewSeq("", nil, alphabet.DNA))
	sc := seqio.NewScanner(r)
	for sc.Next() {
		v = append(v, sc.Seq())
	}
	if sc.Error() != nil {
		log.Fatalf("failed to read sequences: %v", sc.Error())
	}
	//Checking that there aren't too many sequences
	if maxFam != 0 && !subSample && len(v) > maxFam {
		log.Fatalf("too many sequences: %d", len(v))
	}

	var aveLength int
	if subSample {
		// Shuffle first.
		w := make([]seq.Sequence, 0, len(v))
		for _, j := range rand.Perm(len(v)) {
			w = append(w, v[j])
		}
		v = w

		// Calculating lengths and averages for all sequences read
		var LenAllSeq = make([]int, len(v))

		for num, i := range v {
			LenAllSeq[num] = i.Len()
		}
		//	fmt.Printf("The variable: %v\n", LenAllSeq)
		totalLength := 0
		for i := 0; i < len(v); i++ {

			totalLength = totalLength + LenAllSeq[i]
		}
		aveLength = totalLength / len(v)
		fmt.Printf("Average length: %v\n", aveLength)
	}

	// FIXME(brittany): Make this more scientifically robust.
	var (
		sampled int
		buf     bytes.Buffer
	)

	// creating a file for the sequence information
	seqFile := fmt.Sprintf("%s_included_sequences.fq", flag.Args()[0])
	seqOut, sErr := os.Create(filepath.Join(aDir, seqFile))
	if sErr != nil {
		log.Fatalf("failed to create %s: %v", seqFile, sErr)
	}
	defer seqOut.Close()
	seqBufOut := bufio.NewWriter(f)
	defer seqBufOut.Flush()

	// printing subsampled sequences to a file
	var (
		LenSubSeq = make([]int, maxFam) // length of each sequence
	)
	for num, s := range v {
		if sampled++; subSample && sampled > maxFam {
			break
		}
		fmt.Fprintf(&buf, "%60a\n", s)
		fmt.Fprintf(seqOut, "including: %s %s\n", s.Name(), s.Description())

		// Calculating lengths and averages for the subsampled sequences
		LenSubSeq[num] = s.Len()
	}
	totalSubLength := 0
	for i := 0; i < maxFam; i++ {
		totalSubLength = totalSubLength + LenSubSeq[i]
	}
	fmt.Printf("total length of subbed: %v\n", totalSubLength)
	aveSubLength := totalSubLength / maxFam
	fmt.Printf("Average length of subbed: %v\n", aveSubLength)
	fmt.Printf("The sub seq variable for %v sequences: %v\n", maxFam, LenSubSeq)

	var (
		c   *linear.QSeq
		m   *multi.Multi
		err error
	)
	//Creating the consensus

	fmt.Println("Creating consensus")
	c, m, err = consensus(&buf)
	if err != nil {
		log.Printf("failed to generate consensus for %s: %v", flag.Args()[0], err)
		return
	}
	c.ID = fmt.Sprintf("%s_consensus", flag.Args()[0])
	c.Desc = fmt.Sprintf("%d members total %v members sampled", len(v), sampled-1)
	c.Threshold = 42
	c.QFilter = seq.CaseFilter
	conLength := c.Len()

	// Find a way to make the consensus the same case

	// Calculating cutoff of consensus length to mean sequence length
	fmt.Printf("Length of consensus:%v\n", conLength)
	cutoffSubs := float64(conLength) / float64(aveLength)
	cutoffTotal := float64(conLength) / float64(aveSubLength)
	fmt.Printf("ratio for all seqs: %f, \nratio for just sub-sampled: %f\n", cutoffTotal, cutoffSubs)

	// Creating a file for the consensus length information
	confile := fmt.Sprintf("%s_consensus-length-%v.txt", flag.Args()[0], maxFam)
	conOut, err := os.Create(filepath.Join(cDir, confile))
	if err != nil {
		log.Fatalf("failed to create %s: %v", confile, err)
	}

	fmt.Fprintf(conOut, "number sampled: %v \nratio for all seqs: %f \nratio for just sub-sampled: %f\n\n %f\t%f", maxFam, cutoffTotal, cutoffSubs, cutoffTotal, cutoffSubs)

	// Creating a file for the consensus
	file := fmt.Sprintf("%s_consensus.fq%v", flag.Args()[0], maxFam)
	out, err := os.Create(filepath.Join(cDir, file))
	if err != nil {
		log.Fatalf("failed to create %s: %v", file, err)
	}
	if consFasta {
		fmt.Fprintf(out, "%60a\n", c)
	} else {
		fmt.Fprintf(out, "%q\n", c)
	}

	// creating a file for the mutliple alignment
	alignFile := fmt.Sprintf("%s_multiple_alignment%v.fq", flag.Args()[0], maxFam)
	AlignOut, err := os.Create(filepath.Join(aDir, alignFile))
	if err != nil {
		log.Fatalf("failed to create %s: %v", alignFile, err)
	}

	if consFasta {
		fmt.Fprintf(AlignOut, "%60a\n", m)
	} else {
		fmt.Fprintf(AlignOut, "%q\n", m)
	}

	out.Close()

	fmt.Printf("Complete\n\n")
}

func consensus(in io.Reader) (*linear.QSeq, *multi.Multi, error) {
	m, err := muscle.Muscle{Quiet: true}.BuildCommand()
	if err != nil {
		return nil, nil, err
	}
	m.Stdin = in
	buf := &bytes.Buffer{}
	m.Stdout = buf
	err = m.Run()
	if err != nil {
		return nil, nil, err
	}
	var (
		r  = fasta.NewReader(buf, &linear.Seq{Annotation: seq.Annotation{Alpha: alphabet.DNA}})
		ms = &multi.Multi{ColumnConsense: seq.DefaultQConsensus}
	)
	sc := seqio.NewScanner(r)
	for sc.Next() {
		ms.Add(sc.Seq())
	}

	return ms.Consensus(true), ms, sc.Error()
}

func checks(cDir string, aDir string, sDir string) error {

	// check for single fasta
	if len(flag.Args()) != 1 {
		fmt.Fprintln(os.Stderr, "need single input fasta file.")
		flag.Usage()
		os.Exit(1)
	}
	// check, remove and make a consensus target directory
	if cDir == "" {
		fmt.Fprintln(os.Stderr, "need consensus target directory.")
		flag.Usage()
		os.Exit(1)
	}
	cErr := os.RemoveAll(cDir)
	if cErr != nil {
		log.Fatalf("failed to remove consensus target directory: %v", cErr)
		return cErr
	}
	cErr = os.Mkdir(cDir, os.ModeDir|0750)
	if cErr != nil {
		log.Fatalf("failed to create consensus target directory: %v", cErr)
		return cErr
	}

	// check, remove and make an alignment target directory
	if aDir == "" {
		fmt.Fprintln(os.Stderr, "need alignment target directory.")
		flag.Usage()
		os.Exit(1)
	}
	aErr := os.RemoveAll(aDir)
	if aErr != nil {
		log.Fatalf("failed to remove alignment target directory: %v", aErr)
		return aErr
	}
	aErr = os.Mkdir(aDir, os.ModeDir|0750)
	if aErr != nil {
		log.Fatalf("failed to create alignment target directory: %v", aErr)
		return aErr
	}
	// check, remove and make a sequence information target directory

	if sDir == "" {
		fmt.Fprintln(os.Stderr, "need sequence information target directory.")
		flag.Usage()
		os.Exit(1)
	}
	sErr := os.RemoveAll(sDir)
	if sErr != nil {
		log.Fatalf("failed to remove sequence information target directory: %v", sErr)
		return sErr
	}
	sErr = os.Mkdir(sDir, os.ModeDir|0750)
	if sErr != nil {
		log.Fatalf("failed to create sequence information target directory: %v", sErr)
		return sErr
	}

	return nil

}
