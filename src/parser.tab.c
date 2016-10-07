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

void validate();
void realloc_u(int sz);
void realloc_n(int sz);
real get_nearest(const real* invtrans, const real* pos, int* u, int* x, int* y, int* z);
void load_positions(const char* posfilename);
void load_skyrmion(const char* spinsfilename, real* image);
void set_uniform(real* dir, real* spin);
void cut_by_plane(real* o, real* n);
void cut_sphere(real* o, real r);
void append_anisotropy(real scalar, real k[3]);
real* allocate_image();
real apply(const char* name, real arg);
real get_var(const char* name);
int get_var_int(const char* name);

int capacityu, capacityn, dmv_size, ec_size, capacity_anisotropy;

#define YYDEBUG 1
#include "parser.tab.h"

extern FILE *yyin;

#define YYPRINT(file, type, value) yyprint (file, type, value)
static void yyprint(FILE* file, int type, YYSTYPE value); 


#line 104 "src/parser.tab.c" /* yacc.c:339  */

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
    SECTEMP = 278,
    SECCUT = 279,
    PLANE = 280,
    SPHERE = 281,
    VERTEX = 282,
    UNIFORM = 283,
    RANDOM = 284,
    BCFREE = 285,
    BCPERIODIC = 286
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 46 "src/parser.y" /* yacc.c:355  */

 uint sz;
 int i;
 real r;
 real vec[3];
 char* fn;

#line 184 "src/parser.tab.c" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_SRC_PARSER_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 201 "src/parser.tab.c" /* yacc.c:358  */

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
#define YYLAST   198

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  41
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  16
/* YYNRULES -- Number of rules.  */
#define YYNRULES  61
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  154

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   286

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      39,    40,     5,     3,    36,     4,     2,     6,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    37,     2,    38,     2,     2,     2,     2,
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
       2,     2,     2,     2,     2,     2,     1,     2,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint8 yyrline[] =
{
       0,    76,    76,    77,    80,    83,    86,    87,    88,    93,
      98,   103,   109,   110,   111,   112,   113,   113,   117,   123,
     128,   129,   133,   135,   136,   138,   139,   142,   143,   149,
     150,   159,   169,   170,   176,   177,   182,   188,   189,   191,
     194,   195,   199,   204,   209,   210,   216,   217,   218,   220,
     221,   222,   223,   224,   225,   226,   227,   228,   229,   230,
     232,   234
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 1
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "'+'", "'-'", "'*'", "'/'", "NEG",
  "REAL", "SZ", "INTEGER", "FILENAME", "ID", "EOL", "SECS", "SECEF",
  "SECMA", "SECBC", "SECTV", "SECUC", "SECN", "SECEC", "SECDMV",
  "SECIMAGE", "SECLOADIMAGE", "SECPOSITIONS", "SECDIPOLE", "SECTEMP",
  "SECCUT", "PLANE", "SPHERE", "VERTEX", "UNIFORM", "RANDOM", "BCFREE",
  "BCPERIODIC", "','", "'{'", "'}'", "'('", "')'", "$accept",
  "latticefile", "section", "$@1", "bc", "malist", "uclist", "nlist",
  "eclist", "dmvlist", "cutlist", "imagelist", "sz", "integer", "exp",
  "vector", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,    43,    45,    42,    47,   258,   259,   260,
     261,   262,   263,   264,   265,   266,   267,   268,   269,   270,
     271,   272,   273,   274,   275,   276,   277,   278,   279,   280,
     281,   282,   283,   284,   285,   286,    44,   123,   125,    40,
      41
};
# endif

#define YYPACT_NINF -104

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-104)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
    -104,   139,  -104,     5,    27,    43,    47,    53,    70,    71,
      72,    77,    78,    82,    88,    89,    91,    98,  -104,    49,
      52,  -104,   -18,    52,  -104,  -104,  -104,  -104,  -104,   121,
     123,    66,    66,  -104,  -104,  -104,    99,    66,   125,    66,
    -104,  -104,   -18,   132,    52,   109,    66,    55,  -104,   134,
     136,    66,  -104,  -104,  -104,    97,    66,   120,   124,    20,
      49,   104,  -104,    94,    -4,    52,   137,    32,   138,    94,
     163,   -29,  -104,  -104,  -104,    66,    22,    66,    66,    66,
      66,  -104,  -104,    52,    52,    -7,    66,   164,  -104,   166,
     167,  -104,  -104,  -104,  -104,   145,  -104,   170,  -104,    52,
      52,   171,    33,  -104,    63,    63,  -104,  -104,    52,    66,
    -104,    49,    76,  -104,  -104,    52,    32,  -104,    66,   172,
    -104,  -104,   173,   165,   174,    66,  -104,   175,   -28,   112,
    -104,  -104,  -104,  -104,    17,  -104,    32,   180,    66,  -104,
     152,   155,   116,   183,   184,    66,   158,   182,   169,   187,
    -104,  -104,   185,  -104
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     3,     0,
       0,    25,     0,     0,    27,    29,    32,    34,    16,     0,
       0,     0,     0,    37,    44,    45,     0,     0,     0,     7,
      23,    24,     0,     0,    12,    13,    14,    15,    40,     0,
       0,     0,    57,    59,    58,    50,     0,     0,     0,    22,
       0,     0,     6,     0,     0,     0,     0,     0,     0,     0,
       0,    17,    18,    19,    55,     0,     0,     0,     0,     0,
       0,    20,    21,     0,     0,     0,     0,     0,     9,     0,
       0,    28,    47,    46,    48,     0,    33,     0,    36,     0,
       0,     0,     0,    56,    51,    52,    53,    54,     0,     0,
       5,     0,     0,    26,     8,    11,     0,    35,     0,     0,
      43,    49,     0,     0,     0,     0,    61,     0,     0,     0,
      42,    38,    39,     4,     0,    10,     0,     0,     0,    60,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      31,    41,     0,    30
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
    -104,  -104,  -104,  -104,   -30,  -104,  -104,  -104,  -104,  -104,
    -104,  -104,   -59,  -103,   -32,   -12
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     1,    18,    48,    42,    39,    44,    45,    46,    47,
      59,    71,    36,    95,    57,    38
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_uint8 yytable[] =
{
      58,    85,    99,   100,   101,    61,   110,    63,   136,    88,
     137,    43,    64,   128,    68,    69,    40,    41,    19,    74,
      77,    78,    79,    80,    76,    77,    78,    79,    80,   111,
      40,    41,    66,   140,    89,    70,    77,    78,    79,    80,
      20,    92,    93,   102,    94,   104,   105,   106,   107,    83,
      84,    87,   124,    90,   112,   139,    21,    97,    34,    51,
      22,    35,   103,    52,    53,    54,    23,    55,    79,    80,
      51,   108,   109,   121,    52,    53,    54,   123,    55,    77,
      78,    79,    80,    24,    25,    26,   129,   118,   119,    37,
      27,    28,    37,   134,    56,    29,   122,    77,    78,    79,
      80,    30,    31,   127,    32,    56,   142,    77,    78,    79,
      80,    33,   125,   148,   126,    77,    78,    79,    80,    77,
      78,    79,    80,    77,    78,    79,    80,    77,    78,    79,
      80,    37,    49,    81,    50,    60,    75,    82,    62,     2,
      86,    77,    78,    79,    80,    65,    67,    72,   138,    73,
      91,    96,   145,     3,     4,     5,     6,     7,     8,     9,
      10,    11,    12,    13,    14,    15,    16,    17,    77,    78,
      79,    80,    77,    78,    79,    80,    98,   113,   132,   114,
     115,   116,   151,   117,   120,   130,   131,   133,   135,   141,
     143,   144,   146,   147,   149,   150,   152,     0,   153
};

static const yytype_int16 yycheck[] =
{
      32,    60,    31,    32,    33,    37,    13,    39,    36,    13,
      38,    23,    42,   116,    46,    47,    34,    35,    13,    51,
       3,     4,     5,     6,    56,     3,     4,     5,     6,    36,
      34,    35,    44,   136,    64,    47,     3,     4,     5,     6,
      13,     9,    10,    75,    12,    77,    78,    79,    80,    29,
      30,    63,   111,    65,    86,    38,    13,    69,     9,     4,
      13,    12,    40,     8,     9,    10,    13,    12,     5,     6,
       4,    83,    84,    40,     8,     9,    10,   109,    12,     3,
       4,     5,     6,    13,    13,    13,   118,    99,   100,    37,
      13,    13,    37,   125,    39,    13,   108,     3,     4,     5,
       6,    13,    13,   115,    13,    39,   138,     3,     4,     5,
       6,    13,    36,   145,    38,     3,     4,     5,     6,     3,
       4,     5,     6,     3,     4,     5,     6,     3,     4,     5,
       6,    37,    11,    13,    11,    36,    39,    13,    13,     0,
      36,     3,     4,     5,     6,    13,    37,    13,    36,    13,
      13,    13,    36,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    27,    28,     3,     4,
       5,     6,     3,     4,     5,     6,    13,    13,    13,    13,
      13,    36,    13,    13,    13,    13,    13,    13,    13,     9,
      38,    36,     9,     9,    36,    13,     9,    -1,    13
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    42,     0,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    43,    13,
      13,    13,    13,    13,    13,    13,    13,    13,    13,    13,
      13,    13,    13,    13,     9,    12,    53,    37,    56,    46,
      34,    35,    45,    56,    47,    48,    49,    50,    44,    11,
      11,     4,     8,     9,    10,    12,    39,    55,    55,    51,
      36,    55,    13,    55,    45,    13,    56,    37,    55,    55,
      56,    52,    13,    13,    55,    39,    55,     3,     4,     5,
       6,    13,    13,    29,    30,    53,    36,    56,    13,    45,
      56,    13,     9,    10,    12,    54,    13,    56,    13,    31,
      32,    33,    55,    40,    55,    55,    55,    55,    56,    56,
      13,    36,    55,    13,    13,    13,    36,    13,    56,    56,
      13,    40,    56,    55,    53,    36,    38,    56,    54,    55,
      13,    13,    13,    13,    55,    13,    36,    38,    36,    38,
      54,     9,    55,    38,    36,    36,     9,     9,    55,    36,
      13,    13,     9,    13
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    41,    42,    42,    43,    43,    43,    43,    43,    43,
      43,    43,    43,    43,    43,    43,    44,    43,    43,    43,
      43,    43,    43,    45,    45,    46,    46,    47,    47,    48,
      48,    48,    49,    49,    50,    50,    50,    51,    51,    51,
      52,    52,    52,    52,    53,    53,    54,    54,    54,    55,
      55,    55,    55,    55,    55,    55,    55,    55,    55,    55,
      56,    56
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     8,     6,     4,     3,     6,     5,
       8,     6,     3,     3,     3,     3,     0,     4,     4,     4,
       4,     4,     3,     1,     1,     0,     4,     0,     3,     0,
      12,    10,     0,     3,     0,     4,     3,     0,     5,     5,
       0,     9,     4,     3,     1,     1,     1,     1,     1,     4,
       1,     3,     3,     3,     3,     2,     3,     1,     1,     1,
       7,     5
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
#line 80 "src/parser.y" /* yacc.c:1646  */
    { 
	  	sizex=(yyvsp[-5].sz); sizey=(yyvsp[-3].sz); sizez=(yyvsp[-1].sz); 
	  	}
#line 1385 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 5:
#line 83 "src/parser.y" /* yacc.c:1646  */
    { 
	  	sizex=(yyvsp[-3].sz); sizey=(yyvsp[-1].sz); sizez=1; 
	  	}
#line 1393 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 6:
#line 86 "src/parser.y" /* yacc.c:1646  */
    { copy3((yyvsp[-1].vec), magnetic_field); }
#line 1399 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 8:
#line 88 "src/parser.y" /* yacc.c:1646  */
    { 
		boundary_conditions[0]=(yyvsp[-3].sz); 
		boundary_conditions[1]=(yyvsp[-2].sz); 
		boundary_conditions[2]=(yyvsp[-1].sz); 
		}
#line 1409 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 9:
#line 93 "src/parser.y" /* yacc.c:1646  */
    { 
		boundary_conditions[0]=(yyvsp[-2].sz); 
		boundary_conditions[1]=(yyvsp[-1].sz); 
		boundary_conditions[2]=0; 
		}
#line 1419 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 10:
#line 98 "src/parser.y" /* yacc.c:1646  */
    { 
		copy3((yyvsp[-5].vec), translation_vectors[0]);
		copy3((yyvsp[-3].vec), translation_vectors[1]);
		copy3((yyvsp[-1].vec), translation_vectors[2]);
		}
#line 1429 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 11:
#line 103 "src/parser.y" /* yacc.c:1646  */
    { 
		copy3((yyvsp[-3].vec), translation_vectors[0]);
		copy3((yyvsp[-1].vec), translation_vectors[1]);
		real vec[3]={0,0,1};
		copy3(vec, translation_vectors[2]);
		}
#line 1440 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 16:
#line 113 "src/parser.y" /* yacc.c:1646  */
    { 
		real* image=allocate_image(); 
		set_to_field(image); normalize(image);
		}
#line 1449 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 18:
#line 117 "src/parser.y" /* yacc.c:1646  */
    {
		if(!positions) yyerror("Atom positions should be loaded befor image\n"); 
		real* image=allocate_image(); 
		load_skyrmion((yyvsp[-1].fn), image);
		free((yyvsp[-1].fn));
		}
#line 1460 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 19:
#line 123 "src/parser.y" /* yacc.c:1646  */
    {
		if(positions) yyerror("Not unique atom positions definition\n"); 
		load_positions((yyvsp[-1].fn));
		free((yyvsp[-1].fn));
		}
#line 1470 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 20:
#line 128 "src/parser.y" /* yacc.c:1646  */
    { dipole=(yyvsp[-1].r); }
#line 1476 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 21:
#line 129 "src/parser.y" /* yacc.c:1646  */
    { 
		if(temperature!=0) yyerror("Temperature is set twice\n"); 
		temperature=(yyvsp[-1].r); 
		}
#line 1485 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 23:
#line 135 "src/parser.y" /* yacc.c:1646  */
    { (yyval.sz)=0; }
#line 1491 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 24:
#line 136 "src/parser.y" /* yacc.c:1646  */
    { (yyval.sz)=1; }
#line 1497 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 26:
#line 139 "src/parser.y" /* yacc.c:1646  */
    { 
		append_anisotropy((yyvsp[-2].r), (yyvsp[-1].vec)); }
#line 1504 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 28:
#line 143 "src/parser.y" /* yacc.c:1646  */
    {
		if(sizeu>=capacityu) { capacityu=capacityu*2+1; realloc_u(capacityu); };
		copy3((yyvsp[-1].vec), atom_positions+3*sizeu);
		sizeu++;
		}
#line 1514 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 30:
#line 150 "src/parser.y" /* yacc.c:1646  */
    {
		if(sizen>=capacityn) { capacityn=capacityn*2+1; realloc_n(capacityn); };
		neighbours[5*sizen+0]=(yyvsp[-9].i);
		neighbours[5*sizen+1]=(yyvsp[-7].i);
		neighbours[5*sizen+2]=(yyvsp[-5].i);
		neighbours[5*sizen+3]=(yyvsp[-3].sz);
		neighbours[5*sizen+4]=(yyvsp[-1].sz);
		sizen++;
		}
#line 1528 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 31:
#line 159 "src/parser.y" /* yacc.c:1646  */
    {
		if(sizen>=capacityn) { capacityn=capacityn*2+1; realloc_n(capacityn); };
		neighbours[5*sizen+0]=(yyvsp[-7].i);
		neighbours[5*sizen+1]=(yyvsp[-5].i);
		neighbours[5*sizen+2]=0;
		neighbours[5*sizen+3]=(yyvsp[-3].sz);
		neighbours[5*sizen+4]=(yyvsp[-1].sz);
		sizen++;
		}
#line 1542 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 33:
#line 170 "src/parser.y" /* yacc.c:1646  */
    {
		if(ec_size>=capacityn) { capacityn=capacityn*2+1; realloc_n(capacityn); };
		exchange_constant[ec_size]=(yyvsp[-1].r);
		ec_size++;
		}
#line 1552 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 35:
#line 177 "src/parser.y" /* yacc.c:1646  */
    {
		if(dmv_size>=capacityn) { capacityn=capacityn*2+1; realloc_n(capacityn); };
		mult3((yyvsp[-2].r),(yyvsp[-1].vec),dzyaloshinskii_moriya_vector+3*dmv_size);
		dmv_size++;
		}
#line 1562 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 36:
#line 182 "src/parser.y" /* yacc.c:1646  */
    {
		if(dmv_size>=capacityn) { capacityn=capacityn*2+1; realloc_n(capacityn); };
		copy3((yyvsp[-1].vec), dzyaloshinskii_moriya_vector+3*dmv_size);
		dmv_size++;
		}
#line 1572 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 38:
#line 189 "src/parser.y" /* yacc.c:1646  */
    { 
		cut_by_plane((yyvsp[-2].vec),(yyvsp[-1].vec)); }
#line 1579 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 39:
#line 191 "src/parser.y" /* yacc.c:1646  */
    { 
		cut_sphere((yyvsp[-2].vec),(yyvsp[-1].r)); }
#line 1586 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 41:
#line 195 "src/parser.y" /* yacc.c:1646  */
    {
		real* image=initial_state+SIZE*3*(initial_states_count-1);
		append_skyrmion((yyvsp[-6].vec), (yyvsp[-5].r), (yyvsp[-3].r), (yyvsp[-1].r), image);
		}
#line 1595 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 42:
#line 199 "src/parser.y" /* yacc.c:1646  */
    {
		real* image=initial_state+SIZE*3*(initial_states_count-1);
		normalize3((yyvsp[-1].vec));
		set_uniform((yyvsp[-1].vec), image);
		}
#line 1605 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 43:
#line 204 "src/parser.y" /* yacc.c:1646  */
    {
		real* image=initial_state+SIZE*3*(initial_states_count-1);
		skyrmion_random(image);
		}
#line 1614 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 45:
#line 210 "src/parser.y" /* yacc.c:1646  */
    { 
		int n=get_var_int((yyvsp[0].fn)); 
		if(n<0) yyerror("Variable " COLOR_RED "'%s'" COLOR_RESET " value " COLOR_RED "%d" COLOR_RESET " is negative\n", (yyvsp[0].fn), n);
		(yyval.sz)=n;
		}
#line 1624 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 46:
#line 216 "src/parser.y" /* yacc.c:1646  */
    { (yyval.i)=(yyvsp[0].i); }
#line 1630 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 47:
#line 217 "src/parser.y" /* yacc.c:1646  */
    { (yyval.i)=(int)(yyvsp[0].sz); }
#line 1636 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 48:
#line 218 "src/parser.y" /* yacc.c:1646  */
    { (yyval.i)=get_var_int((yyvsp[0].fn)); }
#line 1642 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 49:
#line 220 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=apply((yyvsp[-3].fn), (yyvsp[-1].r)); }
#line 1648 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 50:
#line 221 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=get_var((yyvsp[0].fn)); }
#line 1654 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 51:
#line 222 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=(yyvsp[-2].r)+(yyvsp[0].r); }
#line 1660 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 52:
#line 223 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=(yyvsp[-2].r)-(yyvsp[0].r); }
#line 1666 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 53:
#line 224 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=(yyvsp[-2].r)*(yyvsp[0].r); }
#line 1672 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 54:
#line 225 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=(yyvsp[-2].r)/(yyvsp[0].r); }
#line 1678 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 55:
#line 226 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=-(yyvsp[0].r); }
#line 1684 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 56:
#line 227 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=(yyvsp[-1].r); }
#line 1690 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 58:
#line 229 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=(yyvsp[0].i); }
#line 1696 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 59:
#line 230 "src/parser.y" /* yacc.c:1646  */
    { (yyval.r)=(int)(yyvsp[0].sz); }
#line 1702 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 60:
#line 232 "src/parser.y" /* yacc.c:1646  */
    { 
		(yyval.vec)[0]=(yyvsp[-5].r); (yyval.vec)[1]=(yyvsp[-3].r); (yyval.vec)[2]=(yyvsp[-1].r); }
#line 1709 "src/parser.tab.c" /* yacc.c:1646  */
    break;

  case 61:
#line 234 "src/parser.y" /* yacc.c:1646  */
    { 
		(yyval.vec)[0]=(yyvsp[-3].r); (yyval.vec)[1]=(yyvsp[-1].r); (yyval.vec)[2]=0; }
#line 1716 "src/parser.tab.c" /* yacc.c:1646  */
    break;


#line 1720 "src/parser.tab.c" /* yacc.c:1646  */
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
#line 237 "src/parser.y" /* yacc.c:1906  */


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

real* allocate_image() {
	if(!initial_state) {
		initial_state=(real*)malloc(sizeof(real)*SIZE*3);				
		initial_states_count=1;
	} else {
		initial_states_count++;
		initial_state=(real*)realloc(initial_state, sizeof(real)*SIZE*3*initial_states_count);			
	};
	assert(initial_state);
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
	int empty_mask=active==NULL;
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
		if(empty_mask) SETACTIVE(id)
		else if(!ISACTIVE(id)) { fprintf(stderr, "An attempt to set not active spin at '" COLOR_RED "%s line %d" COLOR_RESET "'\n", spinsfilename, line+2); exit(1); };
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

void set_uniform(real* dir, real* spin) {
	int count=0;
	forall(u,x,y,z) {
		int i=INDEX(u,x,y,z);
		if(ISACTIVE(i)) {
			count++;
			i*=3; for3(j) spin[i+j]=dir[j];
		} else {
			i*=3; for3(j) spin[i+j]=0;
		};
	};
	number_of_active=count;
};

void cut_by_plane(real* o, real* n) {
	int first=active==NULL;
	int count=0;
	real d=dot3(o,n);
	forall(u,x,y,z) {
		real vec[3];
		COORDS(u,x,y,z,vec);
		real a=dot3(vec,n)-d;
		int i=INDEX(u,x,y,z);	
		if(first) {
			if(a>=0) { SETACTIVE(i); count++; }
			else { SETPASSIVE(i); };
		} else {
			if(a>=0) { count++; }
			else { SETPASSIVE(i); };
		};
	};
	number_of_active=count;
};

void cut_sphere(real* o, real r) {
	real rsq=r*r;
	int first=active==NULL;
	int count=0;	
	forall(u,x,y,z) {
		real vec[3];
		COORDS(u,x,y,z,vec);
		sub3(vec,o,vec);
		real a=normsq3(vec)-rsq;
		if(r<0) a=-a;
		int i=INDEX(u,x,y,z);	
		if(first) {
			if(a<=0) { SETACTIVE(i); count++; }
			else { SETPASSIVE(i); };
		} else {
			if(a<=0) { count++; }
			else { SETPASSIVE(i); };
		};
	};
	number_of_active=count;	
};
