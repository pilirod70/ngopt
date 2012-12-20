#include <stdio.h>
#include <stdlib.h>

int compare (const void * a, const void * b) {
  return ( *(int*)a - *(int*)b );
}

void printStats(char* file) {
	FILE * fa_in = fopen (file,"rt");

	char c;
	if (fa_in != NULL) {
		c = (char) getc(fa_in);
		if (c != '>') {
			printf ("File %s not in fasta format\n", file);
			return;
		}
	} else {
		printf ("%s does not exist.\n", file);
		return;
	}
	unsigned int tot = 0;
	unsigned int numN = 0;
	unsigned int numNuc = 0;
	unsigned int len = 0;
	unsigned int numGaps = 0;
	int buf_size = 1024;
	int iter_size = 1024;
	int count = 0;
 	char buf[buf_size];	
	int *tmp_dat;
	unsigned int *dat = malloc(1000*sizeof(int)); 
	int dat_size = 1000;
	int inGap = 0;
	int inHdr = 0;
	int i,j;
	while (!feof(fa_in)){
		iter_size = fread(buf,1,buf_size,fa_in);
		for (i = 0; i < iter_size; i++){
			if (buf[i] == '>'){
				inHdr = 1;
				if (len > 0){
					if (count == dat_size){
						tmp_dat = malloc(2*dat_size*sizeof(int));
						for (j = 0; j < dat_size; j++)
							tmp_dat[j] = dat[j];
						free(dat);
						dat = tmp_dat;
						dat_size = 2*dat_size;
					}
					dat[count++] = len;
				}
				len = 0;
			} else if (buf[i] == '\n') {
				if (inHdr)
					inHdr = 0;
				continue;
			}  else if (!inHdr) {
				if (buf[i] == 'n' || buf[i] == 'N') {
					numN++;
					inGap = 1;
				} else {
					numNuc++;
					if (inGap) {
						inGap = 0;
						numGaps++;	
					}
				}
				len++;
			}
		}
	}

	if (count == dat_size){
		tmp_dat = malloc((dat_size+1)*sizeof(int));
		for (j = 0; j < dat_size; j++)
			tmp_dat[j] = dat[j];
		free(dat);
		dat = tmp_dat;
	}
	dat[count++] = len;

	tot = numN + numNuc;
	qsort(dat, count, sizeof(int), compare);
	unsigned int tally = 0;
	unsigned int n50 = 0;
	unsigned int half = tot/2;
	for (i = 0; i < count; i++){
		tally += dat[i];
		if (tally > half) {
			n50 = dat[i];
			break;
		}
	}
	
	// Fasta_File	Num_Ctgs	N50	AvgCtgLen	MaxCtgLen	MinCtgLen	Total_Bases	Num_Nucs	Num_Unks	NumGaps
	printf("%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",file, count, n50, tot/count, dat[count-1], dat[0], tot, numNuc, numN, numGaps);
}

void usage() {
	printf("Usage: sumfasta [options] <fasta1> <fasta2> ... <fastaN>\n");
	printf("   where options are:\n");
	printf("   --header           print a header indicating the fields of the output\n");
	printf("   -h, --help         print this message and exit\n");
}

int main(int argc, char** argv) {
	if (argc == 1 || (argc == 2 && (strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0))) {
		usage();
		return 1;	
	}
	int start = 1;
	if (strcmp(argv[1],"--header") == 0 ){
		printf("Fasta_File\tNum_Ctgs\tN50\tAvgCtgLen\tMaxCtgLen\tMinCtgLen\tTotal_Bases\tNum_Nucs\tNum_Unks\tNumGaps\n"); 
		start++;
	}
	int i;
	for (i = start; i < argc; i++){
		printStats(argv[i]);
		
	}

}	
