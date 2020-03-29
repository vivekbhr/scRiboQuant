#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <string.h>
//
#include "trimFastq.h"

// THIS PROGRAM WILL TAKE TWO FASTA FILES (R1 & R2) AND TRIM THE UMI sequence OFF R1, adding it in R2

//copy the barcode/UMI sequence (position 0 to 9)
void copy_seq(char d[], char s[]) {
   int bcstrt = 0;
   int bcend = 10;
   //copy first 10 bases
   while (bcstrt < bcend) {
      d[bcstrt] = s[bcstrt];
      bcstrt++;
   }
   d[bcstrt] = '\0';
}

//trim the seq off the barcode and index (10 bases removed)
void trim_seq(char d[], char s[]) {
   int idx = 10;
   int end = strlen(s);
   int i = 0;
   //copy first five bases
   while (idx < end) {
      d[i] = s[idx];
      i++;
      idx++;
   }
   d[i] = '\0';
}

// make header in UMItools format
//void make_header(char * output, char seq[], char tseq[]) {
//  char um[10];
//  char bc[10];
//  output[0] = '\0';
//  copy_seq(um,seq);
//  copy_seq(bc,tseq);
//  strcat(output,bc);
//  strcat(output,"_");
//  strcat(output,um);
//}

// put them together
KSEQ_INIT(gzFile, gzread);

int main(int argc, char *argv[]) //main
{
	gzFile fp;
  gzFile fpt;
	gzFile fpout;
  gzFile fptout;

	kseq_t *seq;
  kseq_t *seqt;

	int l;
	if (argc == 1) {
		fprintf(stderr, "Usage: %s <in_R1.fastq.gz> <in_R2.fastq.gz> <out_R1.fastq> <out_R2.fastq>\n", argv[0]);
		return 1;
	}

	fp = gzopen(argv[1], "r"); // file pointer for R1
  fpt = gzopen(argv[2], "r"); // file pointer for R2

	fpout = gzopen(argv[3], "wb"); // output file pointer for R1
  fptout = gzopen(argv[4], "wb"); // output file pointer for R2

	seq = kseq_init(fp);
  seqt = kseq_init(fpt);
// save data for R1
	char * sequence;
  char * name;
  char * qual;
  char * comment;
// for R2
  char * tsequence;
  char * tname;
  char * tqual;
  char * tcomment;


  char umiseq[10];
  char umiqual[10];
  char barcode[20];
  char trimseq[100];
  char trimqual[100];

	while ((l = kseq_read(seq)) >= 0) {
    kseq_read(seqt);
		sequence = seq->seq.s;
    // other components
    name = seq->name.s;
    if (seq->qual.l) qual = seq->qual.s;//quality
    if (seq->comment.l) comment = seq->comment.s;//comment
    // for file R2
    tsequence = seqt->seq.s;
    tname = seqt->name.s;
    if (seqt->qual.l) tqual = seqt->qual.s;//quality
    if (seqt->comment.l) tcomment = seqt->comment.s;//comment

    //now print data for both R1 and R2

    copy_seq(umiseq,sequence); // copy UMI from R1
    copy_seq(barcode,tsequence); // copy BC from R2
    trim_seq(trimseq,sequence);

    // header = read name + _BC_UMI + comment
    gzprintf(fpout,"@%s %s\n", name, comment);//_%s_%s barcode, umiseq,
    gzprintf(fptout,"@%s %s\n", tname, tcomment);//_%s_%s barcode, umiseq, 

    gzprintf(fpout,"%s\n", trimseq);//trimmed seq for R1
    gzprintf(fptout,"%s%s\n", tsequence, umiseq);//normal seq for R2 + UMI from R1

    gzprintf(fpout,"%s\n","+" );
    gzprintf(fptout,"%s\n","+" );

    copy_seq(umiqual,qual); // copy UMI quality from R1
    trim_seq(trimqual,qual);
    gzprintf(fpout,"%s\n", trimqual);//trimmed qual for R1
    gzprintf(fptout,"%s%s\n", tqual, umiqual);//normal qual for R2 + UMI qual:
	}

	kseq_destroy(seq);
  kseq_destroy(seqt);

	gzclose(fp);
  gzclose(fpt);

  gzclose(fpout);
  gzclose(fptout);

  return 0;

}
