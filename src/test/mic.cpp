
void ffmollac_double(int);

void ffmollac_float(int);

#ifdef GMX_DOUBLE
#define ffmollac ffmollac_double
#else 
#define ffmollac ffmollac_float
#endif // FLOAT

int main()
{
    int a = 0;
    ffmollac(a);
}
