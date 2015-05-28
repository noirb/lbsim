#include <regex.h>
#include <libgen.h>
#include "initLB.h"

int readParameters(
                    int     *xlength,
                    int     *ylength,
                    int     *zlength,
                    double  *tau,
                    double  *velocityWall,
                    int     *timesteps,
                    int     *timestepsPerPlotting,
                    char    *cellDataPath,
                    int     argc,
                    char    *argv[]
                    )
{
  if (argc != 2)
  {
    fprintf(stderr, "ERROR: Incorrect number of inputs given!\n\tUsage:\n\t %s inputFile", argv[0]);
    return -1;
  }

  double velocityWallX, velocityWallY, velocityWallZ;
  char cellDataFile[20];
  char temp[40];

  READ_INT(argv[1],     *xlength);
  READ_INT(argv[1],     *ylength);
  READ_INT(argv[1],     *zlength);
  READ_DOUBLE(argv[1],  *tau);
  READ_DOUBLE(argv[1],  velocityWallX);
  READ_DOUBLE(argv[1],  velocityWallY);
  READ_DOUBLE(argv[1],  velocityWallZ);
  READ_INT(argv[1],     *timesteps);
  READ_INT(argv[1],     *timestepsPerPlotting);
  READ_STRING(argv[1],  cellDataFile);

  velocityWall[0] = velocityWallX;
  velocityWall[1] = velocityWallY;
  velocityWall[2] = velocityWallZ;

  strcpy(temp, argv[1]); // make temporary copy of file path in case dirname modifies the source...
  char *inputFileDir = dirname(temp);
  sprintf(cellDataPath, "%s/%s", inputFileDir, cellDataFile);

  return 0;
}


void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength, int ylength, int zlength, char* cellDataFile)
{
  FILE* cellData;
  char  line[80];
  char  temp[80];
  regex_t regex_single;
  regmatch_t pmatch[10];
  vary_flags wildCardFlags = VARY_NONE;

  const char *pattern_singlecell = "\\((\\*|[0-9]+) (\\*|[0-9]+) (\\*|[0-9]+)\\) ([A-Z]+_*[A-Z]*) *(\\(([0-9]+\\.[0-9]+|[0-9]+) ([0-9]+\\.[0-9]+|[0-9]+) ([0-9]+\\.[0-9]+|[0-9]+)\\))?";
  int reti;

  if( (cellData = fopen(cellDataFile, "r")) == NULL)
  {
    fprintf(stderr, "ERROR: Could not open cell data file: %s\n", cellDataFile);
    exit(1);
  }

        /* Initialize Arrays */
  for (int i = 0; i < xlength+2; i++)
  {
    for (int j = 0; j < ylength+2; j++)
    {
      for (int k = 0; k < zlength+2; k++)
      {
        // set distributions at (i,j,k)
        for (int l = 0; l < NUMBER_OF_LATTICE_DIRECTIONS; l++)
        {
            collideField[INDEXOF(xlength, i, j, k, l)] = LATTICEWEIGHTS[l];
            streamField[INDEXOF(xlength, i, j, k, l)] = LATTICEWEIGHTS[l];
        }

        // set flags at (i,j,k) to FLUID while we're here; obstacle & boundaries set below
        flagField[FINDEXOF(xlength, i, j, k)] = FLUID;
      }
    }
  }

  // read cell data file to set boundary & obstacle conditions
  printf("\n\tParsing cell data file: %s\n", cellDataFile);
  regcomp(&regex_single, pattern_singlecell, REG_EXTENDED);

  int cx, cy, cz, cf;
  while (fgets(line, 80, cellData) != NULL)
  {
    wildCardFlags = VARY_NONE;
    // skip commented lines
    if (line[0] == '#')
    {
        continue;
    }

    // check for cell data on current line & parse it if found
    reti = regexec(&regex_single, line, 10, pmatch, 0);
    if( reti == 0 )
    {
        for (int i = 0; pmatch[i].rm_so != -1; i++)
        {
            int len = pmatch[i].rm_eo - pmatch[i].rm_so;
            strncpy(temp, line + pmatch[i].rm_so, len);
            temp[len] = 0; // null-terminate resulting string
            printf("\t\tMatch %d: '%s'\n", i, temp);
            
            switch(i)                                   /// TODO: Don't do this
            {
                case 0:
                    continue;
                    break;
                case 1:
                    if (strcmp(temp, "*") == 0)
                    {
                        wildCardFlags |= VARY_X;
                        cx = 0;
                    }
                    else
                    {
                        sscanf(temp, "%d", &cx);
                    }
                    break;
                case 2:
                    if (strcmp(temp, "*") == 0)
                    {
                        wildCardFlags |= VARY_Y;
                        cy = 0;
                    }
                    else
                    {
                        sscanf(temp, "%d", &cy);
                    }
                    break;
                case 3:
                    if (strcmp(temp, "*") == 0)
                    {
                        wildCardFlags |= VARY_Z;
                        cz = 0;
                    }
                    else
                    {
                        sscanf(temp, "%d", &cz);
                    }
                    break;
                case 4:
                    if (strcmp(temp, "NO_SLIP") == 0)
                    {
                        cf = NO_SLIP;
                    }
                    else
                    {
                        cf = FLUID;
                    }
                default:
                    continue;
            }
        }

        // if wildcards were found, set all cells that match, otherwise just set this cell
        if (wildCardFlags)
            setFlags(flagField, cf, cx, cy, cz, xlength, ylength, zlength, wildCardFlags);
        else
            flagField[FINDEXOF(xlength, cx, cy, cz)] = cf; // set that flag!
    }
  }
  fclose(cellData);
}

void setFlags(int *flagField, cell_flag flag, int xstart, int ystart, int zstart, int xlength, int ylength, int zlength, vary_flags varying)
{
    int i, max_i; i = max_i = xstart;
    int j, max_j; j = max_j = ystart;
    int k, max_k; k = max_k = zstart;

    if (varying & VARY_X)
    {
        max_i = xlength + 1;
    }
    if (varying & VARY_Y)
    {
        max_j = ylength + 1;
    }
    if (varying & VARY_Z)
    {
        max_k = zlength + 1;
    }

    printf("\t\tSetting wild flags: %d\n\t\t\t x: %d -> %d\n\t\t\t y: %d -> %d\n\t\t\t z: %d -> %d\n", varying, i, max_i, j, max_j, k, max_k);

    for (i = xstart; i <= max_i; i++)
    {
        for (j = ystart; j <= max_j; j++)
        {
            for (k = zstart; k <= max_k; k++)
            {
                flagField[FINDEXOF(xlength, i, j, k)] = flag;
            }
        }
    }
}
