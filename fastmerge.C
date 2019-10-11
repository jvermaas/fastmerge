#include <cstdio>
#include "jsplugin.c"
int main (int argc, char *argv[]) {
	int natoms;
	int i;
	int ntotal = 0;
	if (argc < 3) {
		printf("Requires at least two arguments, an output filename, and the input files used to construct it.\n");
		printf("eg: fastmerge output.js input1.js input2.js input3.js\n");
		exit(0);
	}
	char *outputfilename = argv[1];
	int tbonds = 0;
	int tangles = 0;
	int tdihedral = 0;
	int timproper = 0;
	int tcross = 0;
	long boffset = 0;
	long aoffset = 0;
	long doffset = 0;
	long ioffset = 0;
	long coffset = 0;
	long j;
	//Starts from 2 since the outputfile will be in argv[1]
	for (i=2; i < argc; i++) {
		jshandle *js = (jshandle *) open_js_read(argv[i], "js", &natoms);
		if (js == NULL) {
			//File read failed
			printf("Failed to open file %s\n", argv[i]);
			exit(-1);
		}
		read_js_structure((void *) js, NULL, NULL);
		ntotal += natoms;
		tbonds += js->nbonds;
		tangles += js->numangles;
		tdihedral += js->numdihedrals;
		timproper += js->numimpropers;
		tcross += js->numcterms;
		printf("%d atoms in this structure\n", natoms);
		printf("%d bonds, %d angles, %d dihedrals, %d impropers, %d crossterms\n", js->nbonds, js->numangles, js->numdihedrals, js->numimpropers, js->numcterms);
		close_js_read((void *)js);
	}
	printf("%d atoms in total\n", ntotal);
	printf("%d bonds, %d angles, %d dihedrals, %d impropers, %d crossterms\n", tbonds, tangles, tdihedral, timproper, tcross);
	jshandle *outjs = (jshandle *) open_js_write(outputfilename, "js", ntotal);
	if (outjs == NULL) {
		printf("Failed to open file %s for writing\n", argv[i]);
		exit(-1);
	}
	molfile_atom_t *atoms = (molfile_atom_t*) malloc(sizeof(molfile_atom_t) * ntotal);
	molfile_timestep_t *ts = (molfile_timestep_t*) malloc(sizeof(molfile_timestep_t));
	float *coords = (float *) malloc(3L * ntotal * sizeof(float));
	int *bondfrom = (int *) malloc(tbonds * sizeof(int));
	int *bondto = (int *) malloc(tbonds * sizeof(int));
	int *angles = (int *) malloc(3L * tangles * sizeof(int));
	int *dihedrals = (int *) malloc(4L * tdihedral * sizeof(int));
	int *impropers = (int *) malloc(4L * timproper * sizeof(int));
	int *crossterms = (int *) malloc(8L * tcross * sizeof(int));
	int *optflags = (int*) malloc(sizeof(int));
	int offset = 0;
	for (i=2; i < argc; i++) {
		jshandle *js = (jshandle *) open_js_read(argv[i], "js", &natoms);
		if (js == NULL) {
			//File read failed
			printf("Failed to open file %s\n", argv[i]);
		}
		ts->coords = &(coords[3L*offset]);
		js->verbose = 1;
		read_js_structure((void *) js, optflags, &(atoms[offset]));
		//Add the offset to all the bond indexes.
		for (j = 0; j < js->nbonds; j++) {
			js->bondfrom[j] += offset;
			js->bondto[j] += offset;
		}
		//Copy over bonds from js to the real storage device.
		memcpy(&(bondfrom[boffset]), js->bondfrom, sizeof(int) * js->nbonds);
		memcpy(&(bondto[boffset]), js->bondto, sizeof(int) * js->nbonds);
		boffset += js->nbonds;
		//Now comes the angles!
		for (j = 0; j < 3L * js->numangles; j++) {
			js->angles[j] += offset;
		}
		memcpy(&(angles[aoffset]), js->angles, sizeof(int) * 3L * js->numangles);
		aoffset += 3L * js->numangles;
		//Dihedrals
		for (j = 0; j < 4L * js->numdihedrals; j++) {
			js->dihedrals[j] += offset;
		}
		memcpy(&(dihedrals[doffset]), js->dihedrals, sizeof(int) * 4L * js->numdihedrals);
		doffset += 4L * js->numdihedrals;
		//Impropers
		for (j = 0; j < 4L * js->numimpropers; j++) {
			js->impropers[j] += offset;
		}
		memcpy(&(impropers[ioffset]), js->impropers, sizeof(int) * 4L * js->numimpropers);
		ioffset += 4L * js->numimpropers;
		//Crossterms
		for (j = 0; j < 8L * js->numcterms; j++) {
			js->cterms[j] += offset;
		}
		memcpy(&(crossterms[coffset]), js->cterms, sizeof(int) * 8L * js->numcterms);
		coffset += 8L * js->numcterms;
		//Now read the timestep info.
		read_js_timestep((void *) js, natoms, ts);
		for (j = 0; j < 10; j++) {
			printf("Atom %d: %f %f %f\n", j, coords[3L*(offset+j)+0], coords[3L*(offset+j)+1], coords[3L*(offset+j)+2]);
		}
		offset += natoms;
		close_js_read((void *)js);
	}
	ts->coords=coords;
	outjs->verbose = 1;
	outjs->nbonds = tbonds;
	outjs->bondto = bondto;
	outjs->bondfrom = bondfrom;
	outjs->optflags |= JSOPT_BONDS;
	outjs->optflags |= JSOPT_ANGLES;
	if (tcross > 0) {
		outjs->optflags |= JSOPT_CTERMS;
	}
	outjs->numangles = tangles;
	outjs->angles = angles;
	outjs->numdihedrals = tdihedral;
	outjs->dihedrals = dihedrals;
	outjs->numimpropers = timproper;
	outjs->impropers = impropers;
	outjs->numcterms = tcross;
	outjs->cterms = crossterms;
	write_js_structure((void *)outjs, *optflags, atoms);
	write_js_timestep((void *)outjs, ts);
	printf("%d atoms in total\n", ntotal);
	printf("%d bonds, %d angles, %d dihedrals, %d impropers, %d crossterms\n", outjs->nbonds, outjs->numangles, outjs->numdihedrals, outjs->numimpropers, outjs->numcterms);
	close_js_write((void *)outjs);
	//I'm being sloppy here. I should free stuff, but close_js_write does alot of that for me.
}