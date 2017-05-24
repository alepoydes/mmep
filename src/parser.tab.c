/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
#line 2 "src/parser.y" /* yacc.c:339  */

#include "vector.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>	
#include <string.h>	
#include "parser.h"	
#include "debug.h"	
#include "skyrmion.h"	
#include "assert.h"	
#include "plot.h"	
#include "cut.h"

void validate();
void realloc_u(int sz);
void realloc_n(int sz);
real get_nearest(const real* invtrans, const real* pos, int* u, int* x, int* y, int* z);
void load_positions(const char* posfilename);
void load_skyrmion(const char* spinsfilename, real* image);
void load_from_gnuplot(int need_to_relax, const char* spinsfilename);
void set_uniform(real* dir, real* spin);
//void cut_by_plane(real* o, real* n);
//void cut_sphere(real* o, real r);
void append_anisotropy(real scalar, real k[3]);
real* allocate_image(int);
real apply(const char* name, real arg);
real get_var(const char* name);
int get_var_int(const char* name);
void set_uniform_field(const real vec[3]);
void set_domain(primitive_t* pr);

int capacityu, capacityn, dmv_size, ec_size, capacity_anisotropy;

#define YYDEBUG 1
#include "parser.tab.h"

extern FILE *yyin;

#define YYPRINT(file, type, value) yyprint (file, type, value)
static void yyprint(FILE* file, int type, YYSTYPE value); 


#line 109 "src/parser.tab.c" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 1
#endif

/* In a future release of Bison, this section will be replaced
   by #include "parser.tab.h".  */
#ifndef YY_YY_SRC_PARSER_TAB_H_INCLUDED
# define YY_YY_SRC_PARSER_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    NEG = 258,
    REAL = 259,
    SZ = 260,
    INTEGER = 261,
    FILENAME = 262,
    ID = 263,
    EOL = 264,
    SECS = 265,
    SECEF = 266,
    SECMA = 267,
    SECBC = 268,
    SECTV = 269,
    SECUC = 270,
    SECN = 271,
    SECEC = 272,
    SECDMV = 273,
    SECIMAGE = 274,
    SECLOADIMAGE = 275,
    SECPOSITIONS = 276,
    SECDIPOLE = 277,
    SECIMAGEFROMGNUPLOT = 278,
    SECCURRENT = 279,
    SECTEMP = 280,
    SECCUT = 281,
    SECDOMAIN = 282,
    TIP = 283,
    PLANE = 284,
    SPHERE = 285,
    UNIVERSE = 286,
    VERTEX = 287,
    UNIFORM = 288,
    RANDOM = 289,
    RELAX = 290,
    AFM = 291,
    BCFREE = 292,
    BCPERIODIC = 293
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 53 "src/parser.y" /* yacc.c:355  */

 uint sz;
 int i;
 real r;
 real vec[3];
 char* fn;
 primitive_t* pr;

#line 197 "src/parser.tab.c" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_SRC_PARSER_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 214 "src/parser.tab.c" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   294

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  51
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  20
/* YYNRULES -- Number of rules.  */
#define YYNRULES  86
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  212

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   293

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    10,     2,     2,     2,     2,     4,     2,
      49,    50,     7,     5,    46,     6,     2,     8,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    47,     3,    48,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     9,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    86,    86,    87,    90,    93,    96,    97,    98,   103,
     108,   113,   119,   120,   121,   122,   123,   124,   124,   128,
     134,   139,   144,   145,   149,   153,   159,   160,   162,   163,
     165,   166,   170,   173,   177,   178,   181,   182,   188,   189,
     198,   208,   209,   215,   216,   221,   227,   228,   232,   236,
     237,   241,   246,   247,   251,   256,   260,   265,   266,   267,
     268,   269,   270,   271,   272,   273,   275,   280,   281,   282,
     283,   284,   285,   286,   287,   288,   289,   290,   292,   293,
     294,   298,   299,   300,   301,   303,   305
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 1
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "'|'", "'&'", "'+'", "'-'", "'*'", "'/'",
  "NEG", "'!'", "REAL", "SZ", "INTEGER", "FILENAME", "ID", "EOL", "SECS",
  "SECEF", "SECMA", "SECBC", "SECTV", "SECUC", "SECN", "SECEC", "SECDMV",
  "SECIMAGE", "SECLOADIMAGE", "SECPOSITIONS", "SECDIPOLE",
  "SECIMAGEFROMGNUPLOT", "SECCURRENT", "SECTEMP", "SECCUT", "SECDOMAIN",
  "TIP", "PLANE", "SPHERE", "UNIVERSE", "VERTEX", "UNIFORM", "RANDOM",
  "RELAX", "AFM", "BCFREE", "BCPERIODIC", "','", "'{'", "'}'", "'('",
  "')'", "$accept", "latticefile", "section", "$@1", "is_to_relax", "bc",
  "eflist", "malist", "uclist", "nlist", "eclist", "dmvlist",
  "currentlist", "cutlist", "imagelist", "integer", "sz", "exp", "domain",
  "vector", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   124,    38,    43,    45,    42,    47,   258,
      33,   259,   260,   261,   262,   263,   264,   265,   266,   267,
     268,   269,   270,   271,   272,   273,   274,   275,   276,   277,
     278,   279,   280,   281,   282,   283,   284,   285,   286,   287,
     288,   289,   290,   291,   292,   293,    44,   123,   125,    40,
      41
};
# endif

#define YYPACT_NINF -80

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-80)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     -80,   232,   -80,    -6,     7,    17,    31,    36,    41,    49,
      60,    71,    52,    52,    75,    82,    52,    85,    92,    93,
     102,   -80,    33,   -80,   -80,   -23,    73,   -80,   -80,   -80,
     -80,   -80,   109,   110,   115,    99,   114,   -80,    99,   -80,
     199,    33,   -80,   -80,   -80,    33,   237,    61,    53,    99,
     -80,   -80,   -23,    99,   121,    73,    95,    99,    91,   -80,
     118,   123,    99,   -80,   -80,   -80,    94,    99,   200,   130,
      91,   205,    26,   199,    73,    73,   -80,   199,    14,   -80,
      20,    33,    33,    33,    33,    33,    73,    73,   131,   176,
      15,   180,    73,   133,    33,   262,   176,   134,   -80,   135,
     -80,   -80,    99,    30,    99,    99,    99,    99,   -80,   137,
     176,   138,   -80,    73,    73,   -80,    73,    99,     5,   199,
     199,   -80,   -80,    64,    64,   -80,   -80,   -12,    73,   144,
     -80,   147,   -80,   149,    99,   150,   -80,   184,   -80,   152,
     -80,   -27,   -80,    78,   -80,    67,    67,   -80,   -80,   -80,
     154,   -80,    73,    99,   -80,   282,   -80,   142,   -80,   -80,
      33,   155,   -80,   -80,   -80,   116,    73,    33,   -80,    73,
      73,   159,   181,   -80,   -80,   182,   266,   198,   -80,    99,
     -80,   201,   128,    99,   202,   -80,   -80,   -80,   -80,   -80,
     151,   -80,    33,   192,   188,   -80,   -80,   172,   169,    99,
     207,   210,   195,   178,   209,    99,   215,   -80,   278,   212,
     -80,   -80
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    26,    26,     0,     0,    26,     0,     0,     0,
       0,     3,     0,    30,    34,     0,     0,    36,    38,    41,
      43,    27,     0,     0,     0,     0,     0,    46,     0,    49,
       0,     0,    58,    57,    59,     0,    66,     0,     6,     7,
      28,    29,     0,     0,     0,    13,    14,    15,    16,    17,
       0,     0,     0,    75,    77,    76,    68,     0,     0,     0,
      12,     0,    24,     0,     0,     0,    78,     0,     0,    64,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    52,     0,
      21,    73,     0,     0,     0,     0,     0,     0,    22,     0,
       0,     0,    23,     0,     0,    83,     0,     0,     0,     0,
       0,    25,    65,    60,    61,    62,    63,     0,     0,     0,
      31,     0,     9,     0,     0,     0,    37,     0,    42,     0,
      45,    18,    19,     0,    74,    69,    70,    71,    72,    20,
       0,    48,     0,     0,    79,    80,    84,    81,    82,     5,
       0,     0,    32,    35,     8,     0,    11,     0,    44,     0,
       0,     0,     0,    67,    47,     0,     0,     0,    33,     0,
      86,     0,     0,     0,     0,    55,    56,    50,    51,     4,
       0,    10,     0,     0,     0,    54,    85,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    40,     0,     0,
      53,    39
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -80,   -80,   -80,   -80,    40,   -50,   -80,   -80,   -80,   -80,
     -80,   -80,   -80,   -80,   -80,   -40,   -79,   -38,   -70,     3
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     1,    21,    98,    32,    52,    48,    49,    55,    56,
      57,    58,    70,    72,   141,    46,    47,    68,    78,    54
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_uint8 yytable[] =
{
      71,    79,    90,   115,   159,    80,   127,   118,   119,   120,
      22,    89,   169,   170,   171,    91,   172,   119,   120,    95,
      96,    50,    51,    23,   101,    81,    82,    83,    84,   103,
     121,   132,   110,    24,   160,   104,   105,   106,   107,    41,
     133,   123,   124,   125,   126,    42,    43,    25,    44,   157,
     158,    88,    26,    33,   137,   156,    36,    27,    93,    50,
      51,    97,   113,   114,   143,    28,   145,   146,   147,   148,
     122,    83,    84,   111,   106,   107,    29,   116,   117,   155,
     144,   177,    45,   104,   105,   106,   107,    30,    86,   128,
     129,    34,   131,    87,    31,   135,   165,    62,    35,   139,
      53,    37,    63,    64,    65,    62,    66,    85,    38,    39,
      63,    64,    65,   150,    66,   176,   152,   153,    40,   154,
      53,   104,   105,   106,   107,    59,    60,   182,   173,    61,
      69,   161,    99,    81,    82,    83,    84,    92,    53,   100,
      67,   190,    94,   102,   109,   194,   120,   130,    67,   136,
     140,   142,   197,   149,   151,   175,   104,   105,   106,   107,
     162,   202,   179,   163,   180,   164,   166,   208,   168,   181,
     174,   178,   183,   184,   192,   185,   193,    81,    82,    83,
      84,   104,   105,   106,   107,   104,   105,   106,   107,    81,
      82,    83,    84,   104,   105,   106,   107,   186,   187,   196,
     104,   105,   106,   107,   198,   104,   105,   106,   107,    73,
     104,   105,   106,   107,   189,   201,   108,   191,   195,   203,
     200,   112,   204,    53,   206,   207,   134,   209,   211,     0,
     167,     0,     2,     0,   199,    74,    75,    76,     0,     0,
       0,   205,    81,    82,    83,    84,     0,     0,    77,     3,
       4,     5,     6,     7,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,   104,   105,   106,
     107,   104,   105,   106,   107,     0,     0,     0,   138,     0,
       0,     0,   188,   104,   105,   106,   107,   104,   105,   106,
     107,     0,     0,     0,   210
};

static const yytype_int16 yycheck[] =
{
      38,    41,    52,    73,    16,    45,    85,    77,     3,     4,
      16,    49,    39,    40,    41,    53,    43,     3,     4,    57,
      58,    44,    45,    16,    62,     5,     6,     7,     8,    67,
      16,    16,    70,    16,    46,     5,     6,     7,     8,     6,
      90,    81,    82,    83,    84,    12,    13,    16,    15,   119,
     120,    48,    16,    13,    94,    50,    16,    16,    55,    44,
      45,    58,    36,    37,   102,    16,   104,   105,   106,   107,
      50,     7,     8,    70,     7,     8,    16,    74,    75,   117,
      50,   160,    49,     5,     6,     7,     8,    16,    35,    86,
      87,    16,    89,    40,    42,    92,   134,     6,    16,    96,
      47,    16,    11,    12,    13,     6,    15,    46,    16,    16,
      11,    12,    13,   110,    15,   153,   113,   114,    16,   116,
      47,     5,     6,     7,     8,    16,    16,   167,    50,    14,
      16,   128,    14,     5,     6,     7,     8,    16,    47,    16,
      49,   179,    47,    49,    14,   183,     4,    16,    49,    16,
      16,    16,   192,    16,    16,   152,     5,     6,     7,     8,
      16,   199,    46,    16,    48,    16,    16,   205,    16,   166,
      16,    16,   169,   170,    46,    16,    48,     5,     6,     7,
       8,     5,     6,     7,     8,     5,     6,     7,     8,     5,
       6,     7,     8,     5,     6,     7,     8,    16,    16,    48,
       5,     6,     7,     8,    12,     5,     6,     7,     8,    10,
       5,     6,     7,     8,    16,    46,    16,    16,    16,    12,
      48,    16,    12,    47,    46,    16,    46,    12,    16,    -1,
      46,    -1,     0,    -1,    46,    36,    37,    38,    -1,    -1,
      -1,    46,     5,     6,     7,     8,    -1,    -1,    49,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    33,    34,     5,     6,     7,
       8,     5,     6,     7,     8,    -1,    -1,    -1,    16,    -1,
      -1,    -1,    16,     5,     6,     7,     8,     5,     6,     7,
       8,    -1,    -1,    -1,    16
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    52,     0,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    53,    16,    16,    16,    16,    16,    16,    16,    16,
      16,    42,    55,    55,    16,    16,    55,    16,    16,    16,
      16,     6,    12,    13,    15,    49,    66,    67,    57,    58,
      44,    45,    56,    47,    70,    59,    60,    61,    62,    16,
      16,    14,     6,    11,    12,    13,    15,    49,    68,    16,
      63,    68,    64,    10,    36,    37,    38,    49,    69,    66,
      66,     5,     6,     7,     8,    46,    35,    40,    70,    68,
      56,    68,    16,    70,    47,    68,    68,    70,    54,    14,
      16,    68,    49,    68,     5,     6,     7,     8,    16,    14,
      68,    70,    16,    36,    37,    69,    70,    70,    69,     3,
       4,    16,    50,    66,    66,    66,    66,    67,    70,    70,
      16,    70,    16,    56,    46,    70,    16,    66,    16,    70,
      16,    65,    16,    68,    50,    68,    68,    68,    68,    16,
      70,    16,    70,    70,    70,    68,    50,    69,    69,    16,
      46,    70,    16,    16,    16,    68,    16,    46,    16,    39,
      40,    41,    43,    50,    16,    70,    68,    67,    16,    46,
      48,    70,    66,    70,    70,    16,    16,    16,    16,    16,
      68,    16,    46,    48,    68,    16,    48,    66,    12,    46,
      48,    46,    68,    12,    12,    46,    46,    16,    68,    12,
      16,    16
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    51,    52,    52,    53,    53,    53,    53,    53,    53,
      53,    53,    53,    53,    53,    53,    53,    54,    53,    53,
      53,    53,    53,    53,    53,    53,    55,    55,    56,    56,
      57,    57,    57,    57,    58,    58,    59,    59,    60,    60,
      60,    61,    61,    62,    62,    62,    63,    63,    63,    64,
      64,    64,    65,    65,    65,    65,    65,    66,    66,    66,
      66,    66,    66,    66,    66,    66,    67,    68,    68,    68,
      68,    68,    68,    68,    68,    68,    68,    68,    69,    69,
      69,    69,    69,    69,    69,    70,    70
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     8,     6,     3,     3,     6,     5,
       8,     6,     3,     3,     3,     3,     3,     0,     5,     5,
       5,     4,     4,     4,     3,     4,     0,     1,     1,     1,
       0,     3,     4,     5,     0,     4,     0,     3,     0,    12,
      10,     0,     3,     0,     4,     3,     0,     4,     3,     0,
       5,     5,     0,     9,     4,     3,     3,     1,     1,     1,
       3,     3,     3,     3,     2,     3,     1,     4,     1,     3,
       3,     3,     3,     2,     3,     1,     1,     1,     1,     3,
       3,     3,     3,     2,     3,     7,     5
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                                              );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            /* Fall through.  */
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
{
  YYUSE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yystacksize);

        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 4:
#line 90 "src/parser.y" /* yacc.c:1646  */
    { 
	  	sizex=(yyvsp[-5].sz); sizey=(yyvsp[-3].sz); sizez=(yyvsp[-1].sz); 
	  	}
#line 1446 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 5:
#line 93 "src/parser.y" /* yacc.c:1646  */
    { 
	  	sizex=(yyvsp[-3].sz); sizey=(yyvsp[-1].sz); sizez=1; 
	  	}
#line 1454 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 8:
#line 98 "src/parser.y" /* yacc.c:1646  */
    { 
		boundary_conditions[0]=(yyvsp[-3].sz); 
		boundary_conditions[1]=(yyvsp[-2].sz); 
		boundary_conditions[2]=(yyvsp[-1].sz); 
		}
#line 1464 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 9:
#line 103 "src/parser.y" /* yacc.c:1646  */
    { 
		boundary_conditions[0]=(yyvsp[-2].sz); 
		boundary_conditions[1]=(yyvsp[-1].sz); 
		boundary_conditions[2]=0; 
		}
#line 1474 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 10:
#line 108 "src/parser.y" /* yacc.c:1646  */
    { 
		copy3((yyvsp[-5].vec), translation_vectors[0]);
		copy3((yyvsp[-3].vec), translation_vectors[1]);
		copy3((yyvsp[-1].vec), translation_vectors[2]);
		}
#line 1484 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 11:
#line 113 "src/parser.y" /* yacc.c:1646  */
    { 
		copy3((yyvsp[-3].vec), translation_vectors[0]);
		copy3((yyvsp[-1].vec), translation_vectors[1]);
		real vec[3]={0,0,1};
		copy3(vec, translation_vectors[2]);
		}
#line 1495 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 17:
#line 124 "src/parser.y" /* yacc.c:1646  */
    { 
		real* image=allocate_image((yyvsp[-1].i)); 
		set_to_field(image); normalize(image);
		}
#line 1504 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 19:
#line 128 "src/parser.y" /* yacc.c:1646  */
    {
		if(!positions) yyerror("Atom positions should be loaded befor image\n"); 
		real* image=allocate_image((yyvsp[-3].i)); 
		load_skyrmion((yyvsp[-1].fn), image);
		free((yyvsp[-1].fn));
		}
#line 1515 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 20:
#line 134 "src/parser.y" /* yacc.c:1646  */
    {
		//if(!positions) yyerror("Atom positions should be loaded befor image\n"); 
		load_from_gnuplot((yyvsp[-3].i), (yyvsp[-1].fn));
		free((yyvsp[-1].fn));		
		}
#line 1525 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 21:
#line 139 "src/parser.y" /* yacc.c:1646  */
    {
		if(positions) yyerror("Not unique atom positions definition\n"); 
		load_positions((yyvsp[-1].fn));
		free((yyvsp[-1].fn));
		}
#line 1535 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 22:
#line 144 "src/parser.y" /* yacc.c:1646  */
    { dipole=(yyvsp[-1].r); }
#line 1541 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 23:
#line 145 "src/parser.y" /* yacc.c:1646  */
    { 
		if(temperature!=0) yyerror("Temperature is set twice\n"); 
		temperature=(yyvsp[-1].r); 
		}
#line 1550 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 24:
#line 149 "src/parser.y" /* yacc.c:1646  */
    {
		fprintf(stderr, COLOR_RED"Warning"COLOR_RESET": Section [cut] at %d is "COLOR_RED"depricated"COLOR_RESET". Use "COLOR_BOLD"[domain]"COLOR_RESET" instead\n", yylineno);
		set_domain((yyvsp[0].pr));
		}
#line 1559 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 25:
#line 153 "src/parser.y" /* yacc.c:1646  */
    {
		if(all_active!=NULL) 
			fprintf(stderr, COLOR_RED"Warning"COLOR_RESET" at line %d: domain have been already set.", yylineno);
		set_domain((yyvsp[-1].pr));
		}
#line 1569 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 26:
#line 159 "src/parser.y" /* yacc.c:1646  */
    { (yyval.i)=0; }
#line 1575 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 27:
#line 160 "src/parser.y" /* yacc.c:1646  */
    { (yyval.i)=1; }
#line 1581 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 28:
#line 162 "src/parser.y" /* yacc.c:1646  */
    { (yyval.sz)=0; }
#line 1587 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 29:
#line 163 "src/parser.y" /* yacc.c:1646  */
    { (yyval.sz)=1; }
#line 1593 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 31:
#line 166 "src/parser.y" /* yacc.c:1646  */
    { 
		fprintf(stderr, COLOR_RED "Warning" COLOR_RESET ": Definition of uniform magnetic field should be of the form:\n  " COLOR_BOLD "uniform" COLOR_RESET " {exp,exp,exp}\n");
		set_uniform_field((yyvsp[-1].vec));
		}
#line 1602 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 32:
#line 170 "src/parser.y" /* yacc.c:1646  */
    { 
		set_uniform_field((yyvsp[-1].vec));
		}
#line 1610 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 33:
#line 173 "src/parser.y" /* yacc.c:1646  */
    {
		set_tip_field((yyvsp[-2].vec), (yyvsp[-1].vec));
		}
#line 1618 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 35:
#line 178 "src/parser.y" /* yacc.c:1646  */
    { 
		append_anisotropy((yyvsp[-2].r), (yyvsp[-1].vec)); }
#line 1625 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 37:
#line 182 "src/parser.y" /* yacc.c:1646  */
    {
		if(sizeu>=capacityu) { capacityu=capacityu*2+1; realloc_u(capacityu); };
		copy3((yyvsp[-1].vec), atom_positions+3*sizeu);
		sizeu++;
		}
#line 1635 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 39:
#line 189 "src/parser.y" /* yacc.c:1646  */
    {
		if(sizen>=capacityn) { capacityn=capacityn*2+1; realloc_n(capacityn); };
		neighbours[5*sizen+0]=(yyvsp[-9].i);
		neighbours[5*sizen+1]=(yyvsp[-7].i);
		neighbours[5*sizen+2]=(yyvsp[-5].i);
		neighbours[5*sizen+3]=(yyvsp[-3].sz);
		neighbours[5*sizen+4]=(yyvsp[-1].sz);
		sizen++;
		}
#line 1649 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 40:
#line 198 "src/parser.y" /* yacc.c:1646  */
    {
		if(sizen>=capacityn) { capacityn=capacityn*2+1; realloc_n(capacityn); };
		neighbours[5*sizen+0]=(yyvsp[-7].i);
		neighbours[5*sizen+1]=(yyvsp[-5].i);
		neighbours[5*sizen+2]=0;
		neighbours[5*sizen+3]=(yyvsp[-3].sz);
		neighbours[5*sizen+4]=(yyvsp[-1].sz);
		sizen++;
		}
#line 1663 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 42:
#line 209 "src/parser.y" /* yacc.c:1646  */
    {
		if(ec_size>=capacityn) { capacityn=capacityn*2+1; realloc_n(capacityn); };
		exchange_constant[ec_size]=(yyvsp[-1].r);
		ec_size++;
		}
#line 1673 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 44:
#line 216 "src/parser.y" /* yacc.c:1646  */
    {
		if(dmv_size>=capacityn) { capacityn=capacityn*2+1; realloc_n(capacityn); };
		mult3((yyvsp[-2].r),(yyvsp[-1].vec),dzyaloshinskii_moriya_vector+3*dmv_size);
		dmv_size++;
		}
#line 1683 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 45:
#line 221 "src/parser.y" /* yacc.c:1646  */
    {
		if(dmv_size>=capacityn) { capacityn=capacityn*2+1; realloc_n(capacityn); };
		copy3((yyvsp[-1].vec), dzyaloshinskii_moriya_vector+3*dmv_size);
		dmv_size++;
		}
#line 1693 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 47:
#line 228 "src/parser.y" /* yacc.c:1646  */
    {
		mult3((yyvsp[-2].r),(yyvsp[-1].vec),(yyvsp[-1].vec));
		for3(j) spin_polarized_current[j]+=(yyvsp[-1].vec)[j];
		}
#line 1702 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 48:
#line 232 "src/parser.y" /* yacc.c:1646  */
    {
		for3(j) spin_polarized_current[j]+=(yyvsp[-1].vec)[j];
		}
#line 1710 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 49:
#line 236 "src/parser.y" /* yacc.c:1646  */
    { (yyval.pr)=new_universe(); }
#line 1716 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 50:
#line 237 "src/parser.y" /* yacc.c:1646  */
    { 
		(yyval.pr)=new_intersection((yyvsp[-4].pr), new_plane((yyvsp[-2].vec), (yyvsp[-1].vec)));
		//cut_by_plane($3,$4); 
		}
#line 1725 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 51:
#line 241 "src/parser.y" /* yacc.c:1646  */
    { 
		(yyval.pr)=new_intersection((yyvsp[-4].pr), ((yyvsp[-1].r))>0?new_sphere((yyvsp[-2].vec), (yyvsp[-1].r)):new_complement(new_sphere((yyvsp[-2].vec), -((yyvsp[-1].r)))));
		//cut_sphere($3,$4); 
		}
#line 1734 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 53:
#line 247 "src/parser.y" /* yacc.c:1646  */
    {
		real* image=initial_state+SIZE*3*(initial_states_count-1);
		append_skyrmion((yyvsp[-6].vec), (yyvsp[-5].r), (yyvsp[-3].r), (yyvsp[-1].r), image);
		}
#line 1743 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 54:
#line 251 "src/parser.y" /* yacc.c:1646  */
    {
		real* image=initial_state+SIZE*3*(initial_states_count-1);
		normalize3((yyvsp[-1].vec));
		set_uniform((yyvsp[-1].vec), image);
		}
#line 1753 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 55:
#line 256 "src/parser.y" /* yacc.c:1646  */
    {
		real* image=initial_state+SIZE*3*(initial_states_count-1);
		skyrmion_random(image);
		}
#line 1762 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 56:
#line 260 "src/parser.y" /* yacc.c:1646  */
    {
		real* image=initial_state+SIZE*3*(initial_states_count-1);
		skyrmion_afm(image);
	}
#line 1771 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 57:
#line 265 "src/parser.y" /* yacc.c:1646  */
    { (yyval.i)=(yyvsp[0].i); }
#line 1777 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 58:
#line 266 "src/parser.y" /* yacc.c:1646  */
    { (yyval.i)=(int)(yyvsp[0].sz); }
#line 1783 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 59:
#line 267 "src/parser.y" /* yacc.c:1646  */
    { (yyval.i)=get_var_int((yyvsp[0].fn)); }
#line 1789 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 60:
#line 268 "src/parser.y" /* yacc.c:1646  */
    { (yyval.i)=(yyvsp[-2].i)+(yyvsp[0].i); }
#line 1795 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 61:
#line 269 "src/parser.y" /* yacc.c:1646  */
    { (yyval.i)=(yyvsp[-2].i)-(yyvsp[0].i); }
#line 1801 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 62:
#line 270 "src/parser.y" /* yacc.c:1646  */
    { (yyval.i)=(yyvsp[-2].i)*(yyvsp[0].i); }
#line 1807 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 63:
#line 271 "src/parser.y" /* yacc.c:1646  */
    { (yyval.i)=(yyvsp[-2].i)/(yyvsp[0].i); }
#line 1813 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 64:
#line 272 "src/parser.y" /* yacc.c:1646  */
    { (yyval.i)=-(yyvsp[0].i); }
#line 1819 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 65:
#line 273 "src/parser.y" /* yacc.c:1646  */
    { (yyval.i)=(yyvsp[-1].i); }
#line 1825 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 66:
#line 275 "src/parser.y" /* yacc.c:1646  */
    {
		if((yyvsp[0].i)<=0) yyerror("Value " COLOR_RED "%d" COLOR_RESET " is not positive\n", (yyvsp[0].i));
		(yyval.sz)=(yyvsp[0].i);
	}
#line 1834 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 67:
#line 280 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=apply((yyvsp[-3].fn), (yyvsp[-1].r)); }
#line 1840 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 68:
#line 281 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=get_var((yyvsp[0].fn)); }
#line 1846 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 69:
#line 282 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=(yyvsp[-2].r)+(yyvsp[0].r); }
#line 1852 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 70:
#line 283 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=(yyvsp[-2].r)-(yyvsp[0].r); }
#line 1858 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 71:
#line 284 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=(yyvsp[-2].r)*(yyvsp[0].r); }
#line 1864 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 72:
#line 285 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=(yyvsp[-2].r)/(yyvsp[0].r); }
#line 1870 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 73:
#line 286 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=-(yyvsp[0].r); }
#line 1876 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 74:
#line 287 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=(yyvsp[-1].r); }
#line 1882 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 76:
#line 289 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=(yyvsp[0].i); }
#line 1888 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 77:
#line 290 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=(int)(yyvsp[0].sz); }
#line 1894 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 78:
#line 292 "src/parser.y" /* yacc.c:1646  */
    { (yyval.pr)=new_universe(); }
#line 1900 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 79:
#line 293 "src/parser.y" /* yacc.c:1646  */
    { (yyval.pr)=new_plane((yyvsp[-1].vec), (yyvsp[0].vec)); }
#line 1906 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 80:
#line 294 "src/parser.y" /* yacc.c:1646  */
    { 
		if(!((yyvsp[0].r)>0)) yyerror("Radius should be positive");
		(yyval.pr)=new_sphere((yyvsp[-1].vec), (yyvsp[0].r));
		}
#line 1915 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 81:
#line 298 "src/parser.y" /* yacc.c:1646  */
    { (yyval.pr)=new_union((yyvsp[-2].pr), (yyvsp[0].pr)); }
#line 1921 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 82:
#line 299 "src/parser.y" /* yacc.c:1646  */
    { (yyval.pr)=new_intersection((yyvsp[-2].pr), (yyvsp[0].pr)); }
#line 1927 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 83:
#line 300 "src/parser.y" /* yacc.c:1646  */
    { (yyval.pr)=new_complement((yyvsp[0].pr)); }
#line 1933 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 84:
#line 301 "src/parser.y" /* yacc.c:1646  */
    { (yyval.pr)=(yyvsp[-1].pr); }
#line 1939 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 85:
#line 303 "src/parser.y" /* yacc.c:1646  */
    { 
		(yyval.vec)[0]=(yyvsp[-5].r); (yyval.vec)[1]=(yyvsp[-3].r); (yyval.vec)[2]=(yyvsp[-1].r); }
#line 1946 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 86:
#line 305 "src/parser.y" /* yacc.c:1646  */
    { 
		(yyval.vec)[0]=(yyvsp[-3].r); (yyval.vec)[1]=(yyvsp[-1].r); (yyval.vec)[2]=0; }
#line 1953 "src/parser.tab.c" /* yacc.c:1646  */
    break;


#line 1957 "src/parser.tab.c" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 308 "src/parser.y" /* yacc.c:1906  */


static void yyprint(FILE* file, int type, YYSTYPE value) {
  if (type==FILENAME || type==ID)
    fprintf(file, COLOR_BLUE "%s" COLOR_RESET, value.fn);
  else if(type==SZ)
    fprintf(file, COLOR_BLUE "%d" COLOR_RESET, value.sz);
  else if(type==INTEGER)
    fprintf(file, COLOR_BLUE "%d" COLOR_RESET, value.i);
  else if(type==REAL)
    fprintf(file, COLOR_BLUE "%" RF "g" COLOR_RESET, RT(value.r));
}

void yyerror(const char *s, ...) {
 va_list ap;
 va_start(ap, s);
 fprintf(stderr, COLOR_RED "Parse error" COLOR_RESET " at line %d: ", yylineno);
 vfprintf(stderr, s, ap);
 fprintf(stderr, "\n"); 
 exit(1);
};

void parse_lattice(FILE* file) {
	capacityu=capacityn=sizen=sizeu=dmv_size=ec_size=0; 
	capacity_anisotropy=magnetic_anisotropy_count=0;

	yyin=file;
	yyparse();

	if(number_of_active<0) number_of_active=SIZE;
	validate();
};

real apply(const char* name, real arg) {
	if(strcmp(name,"sin")==0) {
		real c,s; rsincos(arg, &s, &c); return s;
	} else if(strcmp(name,"cos")==0) {
		real c,s; rsincos(arg, &s, &c); return c;
	} else if(strcmp(name,"exp")==0) {
		return rexp(arg);
	} else if(strcmp(name,"print")==0) {
		fprintf(stderr,COLOR_BLUE "%" RF "g" COLOR_RESET "\n", RT(arg));
		return arg;
	};
	yyerror("Undefined function " COLOR_RED "'%s'" COLOR_RESET "\n", name);
	exit(1);
};

int get_var_int(const char* name) {
	const char* val=getenv(name);
	if(val) return atoi(val);
	yyerror("Undefined variable " COLOR_RED "'%s'" COLOR_RESET "\n", name);
	exit(1);	
};

real get_var(const char* name) {
	const char* val=getenv(name);
	if(val) {
		return atof(val);
	} else if(strcmp(name,"pi")==0) {
		return R_PI;
	};
	yyerror("Undefined variable " COLOR_RED "'%s'" COLOR_RESET "\n", name);
	exit(1);	
};

void set_uniform_field(const real vec[3]) {
	if(nonuniform_field) {
		yyerror(COLOR_RED "Error" COLOR_RESET ": Uniform field should be given first\n");
		exit(1);
	};
	if(normsq3(magnetic_field)>0) {
		yyerror(COLOR_RED "Warning" COLOR_RESET ": magnetic field reset\n");
	};
	copy3(vec, magnetic_field); 
};

real* allocate_image(int relax_flag) {
	if(!initial_state) {
		assert(!relax_state);
		initial_state=(real*)malloc(sizeof(real)*SIZE*3);
		initial_states_count=1;
		relax_state=(int*)malloc(sizeof(int));
	} else {
		assert(relax_state);
		initial_states_count++;
		initial_state=(real*)realloc(initial_state, sizeof(real)*SIZE*3*initial_states_count);	
		relax_state=(int*)realloc(relax_state, sizeof(int)*initial_states_count);	
	};
	assert(initial_state); assert(relax_state); 
	relax_state[initial_states_count-1]=relax_flag;
	real* image=initial_state+SIZE*3*(initial_states_count-1);
	return image;
};

void validate() {
	// Check if data is consistent
	if(sizex<0 || sizey<0 || sizez<0) {
		fprintf(stderr,"Parse error: lattice size is negative: %dx%dx%dx%d\n",sizeu,sizex,sizey,sizez);
		exit(1);
	};
	if(SIZE==0) {
		fprintf(stderr,"Parse error: Lattice is empty: %dx%dx%dx%d\n",sizeu,sizex,sizey,sizez);
		exit(1);
	};
	//assert(neighbours); 
	//assert(exchange_constant); 
	//assert(dzyaloshinskii_moriya_vector);
	if(dmv_size!=sizen) {
		fprintf(stderr,"Parse error: Number of DM vectors is %d, expected %d\n",dmv_size,sizen);
		exit(1);
	};
	if(ec_size!=sizen) {
		fprintf(stderr,"Parse error: Number of exchange constants %d, expected %d\n",ec_size,sizen);
		exit(1);
	};
	if(isnan(magnetic_field[0]+magnetic_field[1]+magnetic_field[2])) {
		fprintf(stderr,"Parse error: Magnetic field is not set\n");
		exit(1);
	};
	for(int j=0; j<3; j++) if(boundary_conditions[j]<0 || boundary_conditions[j]>1) {
		fprintf(stderr,"Parse error: Bad boundary conditions: %d,%d,%d\n",boundary_conditions[0],boundary_conditions[1],boundary_conditions[2]);
		exit(1);
	};
	for(uint n=0; n<sizen; n++) {
		if(neighbours[5*n+3]<0 || neighbours[5*n+3]>=sizeu) {
			fprintf(stderr,"Parse error: Wrong source %d for link %d\n",neighbours[5*n+3],n);
			exit(1);
		};
		if(neighbours[5*n+4]<0 || neighbours[5*n+4]>=sizeu) {
			fprintf(stderr,"Parse error: Wrong destination %d for link %d\n",neighbours[5*n+4],n);
			exit(1);
		};
	};
	// Freeing some memory
	realloc_u(sizeu); realloc_n(sizen);
};

void append_anisotropy(real scalar, real k[3]) {
	if(magnetic_anisotropy_count>=capacity_anisotropy) { 
		capacity_anisotropy=capacity_anisotropy*2+1; 
		magnetic_anisotropy=(magnetic_anisotropy_type*)
			realloc(magnetic_anisotropy, sizeof(magnetic_anisotropy_type)*capacity_anisotropy); 
	};
	magnetic_anisotropy[magnetic_anisotropy_count].norm=scalar;
	magnetic_anisotropy[magnetic_anisotropy_count].unit[0]=k[0];
	magnetic_anisotropy[magnetic_anisotropy_count].unit[1]=k[1];
	magnetic_anisotropy[magnetic_anisotropy_count].unit[2]=k[2];
	real t=normsq3(magnetic_anisotropy[magnetic_anisotropy_count].unit);
	magnetic_anisotropy[magnetic_anisotropy_count].norm*=t;
	if(magnetic_anisotropy[magnetic_anisotropy_count].norm>0) {
		magnetic_anisotropy[magnetic_anisotropy_count].unit[0]/=t;
		magnetic_anisotropy[magnetic_anisotropy_count].unit[1]/=t;
		magnetic_anisotropy[magnetic_anisotropy_count].unit[2]/=t;
	};
	magnetic_anisotropy_count++;
};

void realloc_u(int sz) {
	atom_positions=(real*)realloc(atom_positions,sizeof(real)*sz*3);
	assert(atom_positions);
};

void realloc_n(int sz) {
	neighbours=(int*)realloc(neighbours,sizeof(int)*sz*5); assert(neighbours);
	exchange_constant=(real*)realloc(exchange_constant,sizeof(real)*sz); 
	assert(exchange_constant);
	dzyaloshinskii_moriya_vector=(real*)realloc(dzyaloshinskii_moriya_vector,sizeof(real)*sz*3);
	assert(dzyaloshinskii_moriya_vector);
};

real get_nearest(const real* invtrans, const real* pos, int* u, int* x, int* y, int* z) {
	real best=INFINITY;
	for(int lu=0; lu<sizeu; lu++) {
		real p[3]; for3(j) p[j]=pos[j]-atom_positions[lu*3+j];
		real loc[3];
		for3(j) { loc[j]=0; for3(k) loc[j]+=invtrans[3*j+k]*p[k]; };
		int coord[3]; for3(j) coord[j]=round(loc[j]);
		if(coord[0]<0 || coord[0]>=sizex || coord[1]<0 || coord[1]>=sizey || coord[2]<0 || coord[2]>=sizez) continue;
		COORDS(lu,coord[0],coord[1],coord[2],loc);
		sub3(pos, loc, loc);
		real dist=normsq3(loc);
		if(dist<best) {
			best=dist; *u=lu; *x=coord[0]; *y=coord[1]; *z=coord[2];
		};
	};
	return best;
};

void load_positions(const char* posfilename) {
	fprintf(stderr, "Loading positions from '%s'\n", posfilename);
	if(positions) yyerror("Atoms positions should be loaded once\n");
	positions=(int*)malloc(sizeof(int)*SIZE); assert(positions);
	FILE* posfile=fopen(posfilename, "r");
	if(!posfile) { fprintf(stderr, "Can not open " COLOR_RED "%s" COLOR_RESET "\n", posfilename); exit(1); };
	int count; 
   	char buf[256];
   	if(!fgets(buf, sizeof(buf), posfile)) { fprintf(stderr, "There is no header in '" COLOR_RED "%s" COLOR_RESET "'\n", posfilename); exit(1); };
	if(sscanf(buf, "%d", &count)!=1) { fprintf(stderr, "Wrong header of '" COLOR_RED "%s" COLOR_RESET "'\n", posfilename); exit(1); };
	if(count!=SIZE) { fprintf(stderr, "Wrong number of spins in '" COLOR_RED "%s" COLOR_RESET "'\n", posfilename); exit(1); };
	// invert translations matrix
	real invtrans[3][3]; invertmatrix3((real*)translation_vectors, (real*)invtrans);
	real tmp[3][3];	matrixmult3((real*)translation_vectors, (real*)invtrans, (real*)tmp); 
	for3(j) for3(k) if(j==k) assert(rabs(tmp[j][k]-1)<1e-6); else assert(rabs(tmp[j][k])<1e-6);
	for(int line=0; line<SIZE; line++) {
       	if(!fgets(buf, sizeof(buf), posfile)) { fprintf(stderr, "Not enough data in '" COLOR_RED "%s" COLOR_RESET "'\n", posfilename); exit(1); };
		long double posd[3]; 
		int l=sscanf(buf, "%Lg %Lg %Lg", posd, posd+1, posd+2);
		if(l<2 || l>3) { fprintf(stderr, "Position has wrong number of coordinates at '" COLOR_RED "%s line %d" COLOR_RESET "'\n", posfilename, line+2); exit(1); };
		if(l<3) posd[2]=0;
		real pos[3]={posd[0],posd[1],posd[2]};
		int x,y,z,u;
		if(get_nearest((real*)invtrans, pos, &u, &x, &y, &z)>0.01) {
			fprintf(stderr, "Position %" RF "g %" RF "g %" RF "g at '" COLOR_RED "%s line %d" COLOR_RESET "' is too far from lattice\n", RT(pos[0]), RT(pos[1]), RT(pos[2]), posfilename, line+2); 
			exit(1);	
		};
		positions[line]=INDEX(u,x,y,z);
	};
	if(!feof(posfile)) fprintf(stderr, COLOR_YELLOW "Warning:" COLOR_RESET "Extra data in '%s'\n",posfilename);
	fclose(posfile); 
};

void load_skyrmion(const char* spinsfilename, real* image) {
	fprintf(stderr, "Loading image from '%s'\n", spinsfilename);
	// Open files
	FILE* spinsfile=fopen(spinsfilename, "r");
	if(!spinsfile) { 
		fprintf(stderr, "Can not open " COLOR_RED "%s" COLOR_RESET "\n", spinsfilename); 
		exit(1); 
	};
	// Read header
   	char buf[256];
   	if(!fgets(buf, sizeof(buf), spinsfile)) { fprintf(stderr, "There is no header in '" COLOR_RED "%s" COLOR_RESET "'\n", spinsfilename); exit(1); };	
	int count; 
	if(sscanf(buf, "%d", &count)!=1) { fprintf(stderr, "Wrong header of '" COLOR_RED "%s" COLOR_RESET "'\n", spinsfilename); exit(1); };
	if(count!=SIZE) { fprintf(stderr, "Wrong number of spins in '" COLOR_RED "%s" COLOR_RESET "'\n", spinsfilename); exit(1); };	
	// Empty buffer for image
	for(int i=0; i<SIZE*3; i++) image[i]=0.;
	// Read line by line
	int empty_mask=all_active==NULL;
	count=0;
	for(int line=0; line<SIZE; line++) {
       	if(!fgets(buf, sizeof(buf), spinsfile)) { 
       		fprintf(stderr, "Not enough data in '" COLOR_RED "%s" COLOR_RESET "'\n", spinsfilename); 
       		exit(1); 
       	};
		long double spind[3];
		int l=sscanf(buf, "%Lg %Lg %Lg", spind, spind+1, spind+2);
		if(l!=3) { fprintf(stderr, "Spin has wrong number of coordinates at '" COLOR_RED "%s line %d" COLOR_RESET "'\n", spinsfilename, line+2); exit(1); };
		real spin[3]={spind[0],spind[1],spind[2]};
		if(normsq3(spin)<=0.5) continue;
		count++;
		normalize3(spin);
		int id=positions[line];
		if(empty_mask) SETACTIVE(all_active, id)
		else if(!ISACTIVE(all_active, id)) { fprintf(stderr, "An attempt to set not active spin at '" COLOR_RED "%s line %d" COLOR_RESET "'\n", spinsfilename, line+2); exit(1); };
		id*=3;
		for3(j) image[id+j]=spin[j];
	};
	// check if files are consistent
	if(!feof(spinsfile)) fprintf(stderr, COLOR_YELLOW "Warning:" COLOR_RESET "Extra data in '%s'\n",spinsfilename); 
	// Check/set total number of spins
	if(empty_mask) number_of_active=count;
	else if(number_of_active!=count) { 
		fprintf(stderr, "Number of spins %d in image '" COLOR_RED "%s" COLOR_RESET "' does not match number of actvie spins %d\n", count, spinsfilename, number_of_active); 
	exit(1); };
	// finalizing
	fclose(spinsfile);
};
void load_from_gnuplot(int need_to_relax, const char* spinsfilename) {
	real* curried_allocate_image() {
		fprintf(stderr, "  Allocating image\n"); 
		return allocate_image(need_to_relax);
	};
	// Open files
	FILE* spinsfile=fopen(spinsfilename, "r");
	if(!spinsfile) { 
		fprintf(stderr, "Can not open " COLOR_RED "%s" COLOR_RESET "\n", spinsfilename); 
		exit(1); 
	};
	fprintf(stderr, "Reading " COLOR_BLUE "%s" COLOR_RESET "\n", spinsfilename); 
	int count=load_path_from_gnuplot(spinsfile, curried_allocate_image);
	fprintf(stderr, "  Number of loaded images: %d\n", count);
};

void set_uniform(real* dir, real* spin) {
	int count=0;
	forall(u,x,y,z) {
		int i=INDEX(u,x,y,z);
		if(ISACTIVE(all_active, i)) {
			count++;
			i*=3; for3(j) spin[i+j]=dir[j];
		} else {
			i*=3; for3(j) spin[i+j]=0;
		};
	};
	number_of_active=count;
};

void set_domain(primitive_t* pr) {
	fprintf(stderr,COLOR_BLUE"Domain: "COLOR_RESET); fprint_primitive(stderr, pr); fprintf(stderr,"\n");
	int first=all_active==NULL;
	int count=0;
	forall(u,x,y,z) {
		real vec[3];
		COORDS(u,x,y,z,vec);
		bool a=does_belong_to_primitive(pr, vec);
		int i=INDEX(u,x,y,z);	
		if(first) {
			if(a) { SETACTIVE(all_active, i); count++; }
			else { SETPASSIVE(all_active, i); };
		} else {
			if(a) { count++; }
			else { SETPASSIVE(all_active, i); };
		};
	};
	number_of_active=count;
};

/*void cut_by_plane(real* o, real* n) {
	int first=all_active==NULL;
	int count=0;
	real d=dot3(o,n);
	forall(u,x,y,z) {
		real vec[3];
		COORDS(u,x,y,z,vec);
		real a=dot3(vec,n)-d;
		int i=INDEX(u,x,y,z);	
		if(first) {
			if(a>=0) { SETACTIVE(all_active, i); count++; }
			else { SETPASSIVE(all_active, i); };
		} else {
			if(a>=0) { count++; }
			else { SETPASSIVE(all_active, i); };
		};
	};
	number_of_active=count;
};

void cut_sphere(real* o, real r) {
	real rsq=r*r;
	int first=all_active==NULL;
	int count=0;	
	forall(u,x,y,z) {
		real vec[3];
		COORDS(u,x,y,z,vec);
		sub3(vec,o,vec);
		real a=normsq3(vec)-rsq;
		if(r<0) a=-a;
		int i=INDEX(u,x,y,z);	
		if(first) {
			if(a<=0) { SETACTIVE(all_active, i); count++; }
			else { SETPASSIVE(all_active, i); };
		} else {
			if(a<=0) { count++; }
			else { SETPASSIVE(all_active, i); };
		};
	};
	number_of_active=count;	
};*/
