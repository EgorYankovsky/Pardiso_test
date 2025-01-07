#include "stdafx.h"

// MklTestHarm.cpp : Defines the entry point for the console application.
//

ofstream logfile;

int Write_Txt_File_Of_Double(const char *fname, double *massiv, int n_of_records, int len_of_record)
{
	FILE *fp;
	int i, j;

	if ((fp = fopen(fname, "w")) == 0)
	{
		printf("Error: Cannot open file \"%s\" for writing.\n", fname);
		return 1;
	}

	printf("writing %s... ", fname);

	for (i = 0; i < n_of_records; i++)
	{
		for (j = 0; j < len_of_record; j++)
			fprintf(fp, "%25.13e\t", massiv[i*len_of_record + j]);
		fprintf(fp, "\n");
	}
	printf("done\n");

	fclose(fp);
	return 0;
}

int Read_Long_From_Txt_File(const char *fname, int *number)
{
	FILE *fp;
	int temp;
	int retcode;

	if ((fp = fopen(fname, "r")) == 0)
	{
		char str[100];
		sprintf(str, "Error: Cannot open file \"%s\" for reading.\n", fname);
		return 1;
	}

	retcode = fscanf(fp, "%ld", &temp);
	if (retcode != 1)
	{
		char str[100];
		sprintf(str, "Error reading file \"%s\".\n", fname);
		fclose(fp);
	}
	*number = temp;

	fclose(fp);
	return 0;
}

int Read_Bin_File_Of_Double(const char *fname, double *massiv, int n_of_records, int len_of_record)
{
	int temp;
	FILE *fp;

	if ((fp = fopen(fname, "r+b")) == 0)
	{
		char str[100];
		sprintf(str, "Error: Cannot open file \"%s\" for reading.\n", fname);
		return 1;
	}

	temp = fread(massiv, sizeof(double)*len_of_record, n_of_records, fp);
	if (temp != n_of_records)
	{
		char str[100];
		sprintf(str, "Error reading file \"%s\". %ld of %ld records was read.\n", fname, temp, n_of_records);
		fclose(fp);
		return 1;
	}

	fclose(fp);

	return 0;
}

int Read_Bin_File_Of_Long(const char *fname, int *massiv, int n_of_records, int len_of_record)
{
	int temp;
	FILE *fp;

	if ((fp = fopen(fname, "r+b")) == 0)
	{
		char str[100];
		sprintf(str, "Cannot open file %s.\n", fname);
		return 1;
	}

	temp = fread(massiv, sizeof(int)*len_of_record, n_of_records, fp);
	if (temp != n_of_records)
	{
		char str[100];
		sprintf(str, "Error reading file \"%s\". %ld of %ld records was read.\n", fname, temp, n_of_records);
		fclose(fp);
		return 1;
	}


	fclose(fp);
	return 0;
}

void FromRSFToCSR_Real_1_Sym(int nb, int *ig, int *sz_ia, int *sz_ja)
{
	*sz_ia = nb + 1;
	*sz_ja = ig[nb] + nb;
}

void FromRSFToCSR_Real_2_Sym(int nb, int *ig, int *jg, double *di, double *gg,
	MKL_INT *ia, MKL_INT *ja, double *a)
{
	int i, j, k;
	vector<MKL_INT> adr;

	adr.resize(nb, 0);

	for (i = 0; i<nb; i++)
	{
		adr[i] += 1;

		for (j = ig[i]; j <= ig[i + 1] - 1; j++)
		{
			k = jg[j];
			adr[k]++;
		}
	}

	// ia
	ia[0] = 0;
	for (i = 0; i<nb; i++)
		ia[i + 1] = ia[i] + adr[i];

	// ja,  a
	for (i = 0; i<ig[nb] + nb; i++)
		a[i] = 0;

	for (i = 0; i<nb; i++)
		adr[i] = ia[i]; 

	for (i = 0; i<nb; i++)
	{
		ja[adr[i]] = i;
		a[adr[i]] = di[i];
		adr[i]++;
	}

	for (i = 0; i<nb; i++)
	{
		for (j = ig[i]; j <= ig[i + 1] - 1; j++)
		{
			k = jg[j];
			ja[adr[k]] = i;

			a[adr[k]] = gg[j];

			adr[k]++;
		}
	}
}


int _tmain(int argc, _TCHAR* argv[])
{
	logfile.open("pardiso64.log");
	if (!logfile) {
		cerr << "Cannot open pardiso64.log" << endl;
		return 1;
	}

	int i;

	int ig_n_1 = 0;
	int sz_ia = 0;
	int sz_ja = 0;

	int nb;
	int *ig = NULL;
	int *jg = NULL;
	double *di = NULL;
	double *ggl = NULL;

	clock_t begin = clock();


	/* PARDISO() args description */


	/*
		1. void *pt[64];
		   	type: pointers array
		   	description: uses for saving staged data and PARDISO's status
		   	usage: should be declared as zeros array size of 64. 
				   Should stay immutable if matrix's type doesn't changes
	*/
	MKL_INT pt[64]; for (int i(0); i < 64; ++i) pt[i] = 0;

	/*
		2. const MKL_INT *maxfct
		   	type: pointer on integer or long long.
		   	description: maximal factorizations amount.
		   	usage: often sets as 1 to solve one SLAE. 
			recommendation: sets value maxfct = 1.
	*/
	MKL_INT maxfct = 1;


	/*
		3. const MKL_INT *mnum
		   	type: pointer on integer.
		   	description: matrix identifier.
		   	usage: to solve several SLAE sets different mnum's values. 
			recommendation: sets value mnum = 1.
	*/
	MKL_INT	mnum = 1;


	/*
		4. const MKL_INT *mtype
		   	type: pointer on integer.
		   	description: matrix type.
		   	usage: sets matrix type:
				*   2 -> real and symmetric positive definite
				*  11 -> symmetric
				*  -2 -> non symmetric
				*  -4 -> complex and non symmetric
				*   6 -> complex and symmetric
				*  -6 -> complex Hermitian matrix
				* ... -> see Intel MKL documentation for more
	*/
	MKL_INT mtype = 2;


	/*
		5. const MKL_INT *phase
		   	type: pointer on integer.
		   	description: solver's working stage.
		   	usage: sets stage to do:
				*  11 -> analysis and factorization
				*  12 -> numerical factorization
				*  13 -> analysis, factorization and forward substitution
				*  22 -> solve SLAE
				*  -1 -> free memory
			recommendation: use 13 to solve first SLAE, then 22 for each SLAE with this matrix.
						   After usage set phase = -1 to deallocate memory.  
	*/
	MKL_INT phase = 13;


	/*
		6. const MKL_INT *n
		   	type: pointer on integer.
		   	description: matrix size.
	*/
	MKL_INT n = 0;


	/*
		7. void *a
		   	type: array of random type.
		   	description: array of matrix values.
		   	usage: contains matrix values in compressed type:
			recommendation: use float or double types.  
	*/
	double *a = NULL;


	/*
		8. const MKL_INT *ia
		   	type: pointer on array of integers.
		   	description: array of matrix's rows.
		   	usage: contains indexes of rows beginning  
	*/
	MKL_INT *ia = NULL;


	/*
		9. const MKL_INT *ja
		   	type: pointer on array of integers.
		   	description: array of matrix's column.
		   	usage: contains indexes of columns beginning at matrix a  
	*/
	MKL_INT *ja = NULL;


	/*
	   10. const MKL_INT *perm
		   	type: pointer on array of integers.
		   	description: array of permutation (optional).
		   	usage: if permutation is unnecessary perm = nullptr 
	*/
	MKL_INT *perm = NULL;


	/*
	   11. const MKL_INT *nrhs
		   	type: pointer on integer.
		   	description: amount of right-side vectors of SLAE.
		   	usage: if SLAE is only one (Ax=b) then nrhs=1 
	*/
	MKL_INT nrhs = 1;


	/*
	   12. const MKL_INT iparm[64]
		   	type: array of integers.
		   	description: array of PARDISO's parameters.
		   	usage: sets different PARDISO's parameters during solving SLAE
			recommendation: see Intel MKL documentation for more.
	*/
	MKL_INT iparm[64];
	iparm[0] = 1; //iparm(2) - iparm(64) are filled with default values.

				  //iparm[0] = 1; //You must supply all values in components iparm(2) - iparm(64).
	for (i = 1; i < 64; i++) iparm[i] = 0;
	iparm[0] = 1;   // Нет вывода сообщений


	/*
	   13. void *msglvl
		   	type: pointer on void.
		   	description: level of logger output information (optional).
		   	usage: sets value to nullptr if unnecessary
			recommendation: often uses as nullptr. 
						   See Intel MKL documentation for more.
	*/
	MKL_INT msglvl = 1;


	/*
	   14. void *b
		   	type: pointer on void.
		   	description: array of right-side.
		   	usage: type of b should be same as type a 
	*/
	double *pr = NULL;


	/*
	   15. void *x
		   	type: pointer on void.
		   	description: array of solution.
		   	usage: type of x should be same as type b and a 
	*/
	double *x = NULL;


	/*
	   16. MKL_INT *error
		   	type: pointer on integer.
		   	description: variable of error.
		   	usage: if something went wrong PARDISO returns error value at this variable 
	*/
	MKL_INT error = 0;


	// nb
	Read_Long_From_Txt_File("kuslau2", &nb);

	// ig
	ig = new int[nb + 1];
	Read_Bin_File_Of_Long("ig", ig, nb + 1, 1);
	for (i = 0; i<nb + 1; i++) ig[i]--;
	ig_n_1 = ig[nb];

	// jg
	jg = new int[ig_n_1];
	Read_Bin_File_Of_Long("jg", jg, ig_n_1, 1);
	for (i = 0; i<ig_n_1; i++) jg[i]--;

	// di
	di = new double[nb];
	Read_Bin_File_Of_Double("di", di, nb, 1);

	// ggl
	ggl = new double[ig_n_1];
	Read_Bin_File_Of_Double("gg", ggl, ig_n_1, 1);

	// pr
	pr = new double[nb];
	Read_Bin_File_Of_Double("pr", pr, nb, 1);

	FromRSFToCSR_Real_1_Sym(nb, ig, &sz_ia, &sz_ja);

	ia = new MKL_INT[sz_ia];
	ja = new MKL_INT[sz_ja];
	a = new double[sz_ja];

	FromRSFToCSR_Real_2_Sym(nb, ig, jg, di, ggl, ia, ja, a);

	for (i = 0; i<sz_ia; i++) ia[i]++;

	for (i = 0; i<sz_ja; i++) ja[i]++;

	if (ig) { delete[] ig; ig = NULL; }
	if (jg) { delete[] jg; jg = NULL; }
	if (di) { delete[] di; di = NULL; }
	if (ggl) { delete[] ggl; ggl = NULL; }


	// PARDISO
	n = nb;
	x = new double[nb];
	perm = new MKL_INT[nb];

	cout << "pardiso start.." << endl << flush;
	PARDISO(
		pt,			// 1
		&maxfct,	// 2
		&mnum,		// 3
		&mtype,		// 4
		&phase,		// 5
		&n,			// 6
		a,			// 7
		ia,			// 8
		ja,			// 9
		perm,		// 10
		&nrhs,		// 11
		iparm,		// 12
		&msglvl,	// 13
		pr,			// 14
		x,			// 15
		&error		// 16
	);

	phase = -1;
	PARDISO(
		pt, 
		&maxfct, 
		&mnum, 
		&mtype, 
		&phase, 
		&n, 
		a, 
		ia, 
		ja, 
		perm, 
		&nrhs, 
		iparm, 
		&msglvl, 
		pr, 
		x, 
		&error
	);

	clock_t time = (clock() - begin) / CLOCKS_PER_SEC;

	ofstream fout;
	fout.open("kit", ios_base::app);
	fout << "n=" << nb << " pardiso time=" << time << " s" << endl;
	fout.close();

	Write_Txt_File_Of_Double("x.txt", x, nb, 1);

	cout << "info=" << error << endl;
	logfile << "info=" << error << endl;

	if (a) { delete[] a;  a = NULL; }
	if (x) { delete[] x;  x = NULL; }
	if (pr) { delete[] pr; pr = NULL; }
	if (ia) { delete[] ia; ia = NULL; }
	if (ja) { delete[] ja; ja = NULL; }
	if (perm) { delete[] perm; perm = NULL; }
	

	logfile.close();
	logfile.clear();

	return 0;
}