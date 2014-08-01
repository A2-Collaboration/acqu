void Clean()
{
  Char_t Buff[256];
  Char_t Name[256];
  Char_t Line1[1024];
  Char_t Line2[1024];
  FILE* In;
  FILE* Out1;
  FILE* Out2;

  In   = fopen("/home/svenschu/4v4a/data/Uncalib.txt", "r");
  Out1 = fopen("/home/svenschu/4v4a/data/Element0.txt", "w");
  Out2 = fopen("/home/svenschu/4v4a/data/TimeWalk.txt", "w");
  while(!feof(In))
  {
    if(fgets(Line1, 1023, In)) fputs(Line1, Out1);
    if(fgets(Line2, 1023, In)) fputs(Line2, Out2);
  }
  fclose(In);
  fclose(Out1);
  fclose(Out2);
}

void Merge()
{
  Char_t Buff[256];
  Char_t Name[256];
  Char_t Line1[1024];
  Char_t Line2[1024];
  FILE* In1;
  FILE* In2;
  FILE* Out;

  In1  = fopen("/home/svenschu/4v4a/data/Element1.txt", "r");
  In2  = fopen("/home/svenschu/4v4a/data/TimeWalk.txt", "r");
  Out  = fopen("/home/svenschu/4v4a/data/Calib.txt", "w");
  while(!feof(In1))
  {
    if(fgets(Line1, 1023, In1)) fputs(Line1, Out);
    if(fgets(Line2, 1023, In2)) fputs(Line2, Out);
  }
  fclose(In1);
  fclose(In2);
  fclose(Out);
}

void Reformat()
{
  Char_t Buff[256];
  Char_t Name[256];
  Double_t Val[14];
  FILE* In;
  FILE* Out;

  In  = fopen("/home/svenschu/4v4a/data/Unsorted.txt", "r");
  Out = fopen("/home/svenschu/4v4a/data/Sorted.txt", "w");
  while(!feof(In))
  {
    fscanf(In, "Element: %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
           &Val[0],  &Val[1],  &Val[2],  &Val[3],  &Val[4],  &Val[5],  &Val[6],  &Val[7],
           &Val[8],  &Val[9],  &Val[10], &Val[11], &Val[12]);
    fprintf(Out, "Element:\n");
  }
  fclose(In);
  fclose(Out);
}

