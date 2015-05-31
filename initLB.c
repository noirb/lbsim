#include <ctype.h>
#include <regex.h>
#include <libgen.h>
#include "initLB.h"

int readParameters(
                    int     *xlength,
                    int     *ylength,
                    int     *zlength,
                    double  *tau,
                    int     *timesteps,
                    int     *timestepsPerPlotting,
                    int     argc,
                    char    *argv[]
                    )
{
  if (argc != 2)
  {
    fprintf(stderr, "ERROR: Incorrect number of inputs given!\n\tUsage:\n\t %s inputFile", argv[0]);
    return -1;
  }

  READ_INT(argv[1],     *xlength);
  READ_INT(argv[1],     *ylength);
  READ_INT(argv[1],     *zlength);
  READ_DOUBLE(argv[1],  *tau);
  READ_INT(argv[1],     *timesteps);
  READ_INT(argv[1],     *timestepsPerPlotting);

  return 0;
}


void initialiseFields(double *collideField, double *streamField, flag_data *flagField, int xlength, int ylength, int zlength, char* cellDataFile)
{
  FILE* cellData;
  char  line[80];
  char  temp[80];
  regex_t regex_single;
  regmatch_t pmatch[10];
  vary_flags wildCardFlags = VARY_NONE;

  const char *pattern_singlecell = "\\((\\*|[nN]|[0-9]+) +(\\*|[nN]|[0-9]+) +(\\*|[nN]|[0-9]+)\\) +([A-Z]+_*[A-Z]*) *(\\((-?[0-9]+\\.[0-9]+|-?[0-9]+) *(-?[0-9]+\\.[0-9]+|-?[0-9]+)? *(-?[0-9]+\\.[0-9]+|-?[0-9]+)?\\))?";
  int reti;
  double cellParms[3] = {0, 0, 0};

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
            collideField[INDEXOF(i, j, k, l)] = LATTICEWEIGHTS[l];
            streamField[INDEXOF(i, j, k, l)] = LATTICEWEIGHTS[l];
        }

        // set flags at (i,j,k) to FLUID while we're here; obstacle & boundaries set below
        flagField[FINDEXOF(i, j, k)].flag = FLUID;
      }
    }
  }

  // read cell data file to set boundary & obstacle conditions
  printf("\n\tParsing cell data from file: %s\n", cellDataFile);
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

            switch(i)
            {
                case 0:
                    continue;
                    break;
                case 1:
                    if (temp[0] == '*')
                    {
                        wildCardFlags |= VARY_X;
                        cx = 0;
                    }
                    else if (toupper((int)temp[0]) == 'N')
                    {
                        cx = xlength+1;
                    }
                    else
                    {
                        sscanf(temp, "%d", &cx);
                    }
                    break;
                case 2:
                    if (temp[0] == '*')
                    {
                        wildCardFlags |= VARY_Y;
                        cy = 0;
                    }
                    else if (toupper((int)temp[0]) == 'N')
                    {
                        cy = ylength+1;
                    }
                    else
                    {
                        sscanf(temp, "%d", &cy);
                    }
                    break;
                case 3:
                    if (temp[0] == '*')
                    {
                        wildCardFlags |= VARY_Z;
                        cz = 0;
                    }
                    else if (toupper((int)temp[0]) == 'N')
                    {
                        cz = zlength+1;
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
                    else if (strcmp(temp, "MOVING_WALL") == 0)
                    {
                        cf = MOVING_WALL;
                    }
                    else if (strcmp(temp, "INFLOW") == 0)
                    {
                        cf = INFLOW;
                    }
                    else if (strcmp(temp, "OUTFLOW") == 0)
                    {
                        cf = OUTFLOW;
                    }
                    else if (strcmp(temp, "FREE_SLIP") == 0)
                    {
                        cf = FREE_SLIP;
                    }
                    else if (strcmp(temp, "PRESSURE_IN") == 0)
                    {
                        cf = PRESSURE_IN;
                    }
                    else
                    {
                        cf = FLUID;
                    }
                case 5: // contains all parameter groups
                    continue;
                    break;
                case 6: // first parameter
                    sscanf(temp, "%lf", cellParms);
                    break;
                case 7: // second parameter
                    sscanf(temp, "%lf", cellParms+1);
                    break;
                case 8: // third parameter
                    sscanf(temp, "%lf", cellParms+2);
                    break;
                default:
                    continue;
            }
        }

        // if wildcards were found, set all cells that match, otherwise just set this cell
        if (wildCardFlags)
        {
            setFlags(flagField, cf, cx, cy, cz, cellParms, xlength, ylength, zlength, wildCardFlags);
        }
        else
        {
            flagField[FINDEXOF(cx, cy, cz)].flag = cf; // set that flag!
            flagField[FINDEXOF(cx, cy, cz)].parms[0] = cellParms[0];
            flagField[FINDEXOF(cx, cy, cz)].parms[1] = cellParms[1];
            flagField[FINDEXOF(cx, cy, cz)].parms[2] = cellParms[2];
        }
    }
  }
  fclose(cellData);
}

void setFlags(flag_data *flagField, cell_flag flag, int xstart, int ystart, int zstart, double* cell_parameters, int xlength, int ylength, int zlength, vary_flags varying)
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

    for (i = xstart; i <= max_i; i++)
    {
        for (j = ystart; j <= max_j; j++)
        {
            for (k = zstart; k <= max_k; k++)
            {
                flagField[FINDEXOF(i, j, k)].flag = flag;
                flagField[FINDEXOF(i, j, k)].parms[0] = cell_parameters[0];
                flagField[FINDEXOF(i, j, k)].parms[1] = cell_parameters[1];
                flagField[FINDEXOF(i, j, k)].parms[2] = cell_parameters[2];
            }
        }
    }
}
