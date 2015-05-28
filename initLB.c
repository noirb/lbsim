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
  regex_t regex;
  regmatch_t pmatch[10];

  const char *pattern = "\\(([0-9]+) ([0-9]+) ([0-9]+)\\) ([A-Z]+_*[A-Z]*) (\\(([0-9]+\\.[0-9]+|[0-9]+) ([0-9]+\\.[0-9]+|[0-9]+) ([0-9]+\\.[0-9]+|[0-9]+)\\))?";
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
  if (regcomp(&regex, pattern, REG_EXTENDED))
  {
    fprintf(stderr, "ERROR: Could not compile regex for file parsing!\n");
    exit(1);
  }

  int cx, cy, cz, cf;
  while (fgets(line, 80, cellData) != NULL)
  {
    // skip commented lines
    if (line[0] == '#')
    {
        continue;
    }

    // check
    printf("Checking line: '%s'\n", line);
    reti = regexec(&regex, line, 10, pmatch, 0);
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
                    sscanf(temp, "%d", &cx);
                    break;
                case 2:
                    sscanf(temp, "%d", &cy);
                    break;
                case 3:
                    sscanf(temp, "%d", &cz);
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
        flagField[FINDEXOF(xlength, cx, cy, cz)] = cf; // set that flag!
    }
    else
    {
        printf("\t\tNo matches on line: %s", line);
    }
  }
  
}

