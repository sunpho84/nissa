/* A Bison parser, made by GNU Bison 3.2.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018 Free Software Foundation, Inc.

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

/* Undocumented macros, especially those whose name start with YY_,
   are private implementation details.  Do not rely on them.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.2"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1


/* Substitute the variable and function names.  */
#define yyparse         parser_parse
#define yylex           parser_lex
#define yyerror         parser_error
#define yydebug         parser_debug
#define yynerrs         parser_nerrs


/* First part of user prologue.  */
#line 9 "../../projects/generate_confs/parser.ypp" /* yacc.c:338  */


#define YYDEBUG 1

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <glob.h>

#include "nissa.hpp"
#include "generate_confs/driver.hpp"
#include "generate_confs/parser.hpp"

  using namespace nissa;

  const int debug_parser=0;
  
  int tokenizer_lex(YYSTYPE *lvalp,YYLTYPE *llocp,void *scanner);
#define parser_lex tokenizer_lex

  void parser_error(YYLTYPE *locp,driver_t *driver,const char *err)
  {crash("exception at line %d columns [%d-%d] %s",locp->first_line,locp->first_column,locp->last_column,err);}

#define scanner driver->scanner
  

#line 102 "generate_confs/parser.cpp" /* yacc.c:338  */
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
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
   by #include "y.tab.h".  */
#ifndef YY_PARSER_GENERATE_CONFS_PARSER_HPP_INCLUDED
# define YY_PARSER_GENERATE_CONFS_PARSER_HPP_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int parser_debug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    TK_INVALID_STRING = 258,
    TK_INVALID_CHAR = 259,
    TK_DOUBLE = 260,
    TK_INT = 261,
    TK_QUANTITY = 262,
    TK_MINUS = 263,
    TK_PLUS = 264,
    TK_TIMES = 265,
    TK_DIV = 266,
    NEG = 267,
    TK_POW = 268,
    TK_GEOMETRY = 269,
    TK_LX = 270,
    TK_LY = 271,
    TK_LZ = 272,
    TK_L = 273,
    TK_T = 274,
    TK_TAG = 275,
    TK_THEORY = 276,
    TK_GAUGE_PARS = 277,
    TK_BETA = 278,
    TK_QUOTED_TEXT = 279,
    TK_EACH = 280,
    TK_MEAS_EACH_NSMOOTH = 281,
    TK_AFTER = 282,
    TK_NONE = 283,
    TK_YES = 284,
    TK_NO = 285,
    TK_RESIDUE = 286,
    TK_META = 287,
    TK_ORDINARY = 288,
    TK_TOPO_POTENTIAL = 289,
    TK_THETA = 290,
    TK_COEFF = 291,
    TK_WIDTH = 292,
    TK_BARR = 293,
    TK_FORCE_OUT = 294,
    TK_WELL_TEMPERING = 295,
    TK_BEND = 296,
    TK_PATH = 297,
    TK_META_PRINT_SCANNED_INPUT = 298,
    TK_META_PRINT_FULL_INPUT = 299,
    TK_META_PRINT = 300,
    TK_QUARK = 301,
    TK_DEGENERACY = 302,
    TK_DISCRETIZ = 303,
    TK_ROOT_STAG = 304,
    TK_OVERLAP = 305,
    TK_ROOT_TM_CLOV = 306,
    TK_KAPPA = 307,
    TK_CSW = 308,
    TK_MASS = 309,
    TK_MASS_OVERLAP = 310,
    TK_RE_POT_CH = 311,
    TK_IM_POT_CH = 312,
    TK_ELEC_CHARGE = 313,
    TK_SMOOTH_METHOD = 314,
    TK_SPACE_OR_TIME = 315,
    TK_SPACE = 316,
    TK_TIME = 317,
    TK_SPACETIME = 318,
    TK_MEAS_CORR = 319,
    TK_CORR_PATH = 320,
    TK_NOISE_TYPE = 321,
    TK_RND_T = 322,
    TK_STOUT = 323,
    TK_NLEVELS = 324,
    TK_RHO = 325,
    TK_COOLING = 326,
    TK_NSTEPS = 327,
    TK_APE = 328,
    TK_ALPHA = 329,
    TK_HYP = 330,
    TK_ALPHAS = 331,
    TK_WFLOW = 332,
    TK_NFLOWS = 333,
    TK_NRECU = 334,
    TK_FLOW_STEP = 335,
    TK_GAUGE_ACTION = 336,
    TK_WILSON = 337,
    TK_TLSYM = 338,
    TK_IWASAKI = 339,
    TK_BKGRD_EM_FIELD = 340,
    TK_B_COMP = 341,
    TK_E_COMP = 342,
    TK_ITHEORY = 343,
    TK_NCOPIES = 344,
    TK_NHITS = 345,
    TK_MEAS_MESON_CORRS = 346,
    TK_OPERATORS = 347,
    TK_MEAS_NUCLEON_CORRS = 348,
    TK_MEAS_PUTPOURRI = 349,
    TK_COMPUTE_SUSC = 350,
    TK_MAX_ORDER = 351,
    TK_MEAS_RENDENS = 352,
    TK_MEAS_ZUMBA = 353,
    TK_MEAS_SPINPOL = 354,
    TK_USE_FERM_CONF_FOR_GLUONS = 355,
    TK_USE_ADJOINT_FLOW = 356,
    TK_MEAS_QED_CORRS = 357,
    TK_MEAS_MAGNETIZ = 358,
    TK_MEAS_MIN_MAX_EIGENVAL = 359,
    TK_MIN_MAX = 360,
    TK_SPECTR_PROJ = 361,
    TK_NEIGS = 362,
    TK_EIG_PRECISION = 363,
    TK_WSPACE_SIZE = 364,
    TK_USE_SMOOTH = 365,
    TK_MEAS_PLAQ_POL = 366,
    TK_MEAS_PLAQ = 367,
    TK_MEAS_ENERGY = 368,
    TK_MEAS_POLY = 369,
    TK_MEAS_TOP = 370,
    TK_MEAS_LUPPOLI = 371,
    TK_MEAS_WATUSSO = 372,
    TK_MEAS_ALL_RECTS = 373,
    TK_SPATIAL = 374,
    TK_TEMPORAL = 375,
    TK_DMIN = 376,
    TK_DMAX = 377,
    TK_TMIN = 378,
    TK_TMAX = 379,
    TK_EVOLUTION = 380,
    TK_ID_SEA_THEORY = 381,
    TK_NTRAJ_TOT = 382,
    TK_SKIP_METRO = 383,
    TK_TRAJ_LENGTH = 384,
    TK_ACT_RESIDUE = 385,
    TK_MD_RESIDUE = 386,
    TK_NSUBSTEPS = 387,
    TK_NPSEUDO_FERMS = 388,
    TK_GAUGE_CONF = 389,
    TK_STORE_PATH = 390,
    TK_STORE_EACH = 391,
    TK_STORE_RUNNING = 392,
    TK_START_COND = 393,
    TK_HOT = 394,
    TK_COLD = 395,
    TK_ANALYSIS = 396,
    TK_CONF_LIST = 397,
    TK_RUN = 398,
    TK_WALLTIME = 399,
    TK_SEED = 400
  };
#endif
/* Tokens.  */
#define TK_INVALID_STRING 258
#define TK_INVALID_CHAR 259
#define TK_DOUBLE 260
#define TK_INT 261
#define TK_QUANTITY 262
#define TK_MINUS 263
#define TK_PLUS 264
#define TK_TIMES 265
#define TK_DIV 266
#define NEG 267
#define TK_POW 268
#define TK_GEOMETRY 269
#define TK_LX 270
#define TK_LY 271
#define TK_LZ 272
#define TK_L 273
#define TK_T 274
#define TK_TAG 275
#define TK_THEORY 276
#define TK_GAUGE_PARS 277
#define TK_BETA 278
#define TK_QUOTED_TEXT 279
#define TK_EACH 280
#define TK_MEAS_EACH_NSMOOTH 281
#define TK_AFTER 282
#define TK_NONE 283
#define TK_YES 284
#define TK_NO 285
#define TK_RESIDUE 286
#define TK_META 287
#define TK_ORDINARY 288
#define TK_TOPO_POTENTIAL 289
#define TK_THETA 290
#define TK_COEFF 291
#define TK_WIDTH 292
#define TK_BARR 293
#define TK_FORCE_OUT 294
#define TK_WELL_TEMPERING 295
#define TK_BEND 296
#define TK_PATH 297
#define TK_META_PRINT_SCANNED_INPUT 298
#define TK_META_PRINT_FULL_INPUT 299
#define TK_META_PRINT 300
#define TK_QUARK 301
#define TK_DEGENERACY 302
#define TK_DISCRETIZ 303
#define TK_ROOT_STAG 304
#define TK_OVERLAP 305
#define TK_ROOT_TM_CLOV 306
#define TK_KAPPA 307
#define TK_CSW 308
#define TK_MASS 309
#define TK_MASS_OVERLAP 310
#define TK_RE_POT_CH 311
#define TK_IM_POT_CH 312
#define TK_ELEC_CHARGE 313
#define TK_SMOOTH_METHOD 314
#define TK_SPACE_OR_TIME 315
#define TK_SPACE 316
#define TK_TIME 317
#define TK_SPACETIME 318
#define TK_MEAS_CORR 319
#define TK_CORR_PATH 320
#define TK_NOISE_TYPE 321
#define TK_RND_T 322
#define TK_STOUT 323
#define TK_NLEVELS 324
#define TK_RHO 325
#define TK_COOLING 326
#define TK_NSTEPS 327
#define TK_APE 328
#define TK_ALPHA 329
#define TK_HYP 330
#define TK_ALPHAS 331
#define TK_WFLOW 332
#define TK_NFLOWS 333
#define TK_NRECU 334
#define TK_FLOW_STEP 335
#define TK_GAUGE_ACTION 336
#define TK_WILSON 337
#define TK_TLSYM 338
#define TK_IWASAKI 339
#define TK_BKGRD_EM_FIELD 340
#define TK_B_COMP 341
#define TK_E_COMP 342
#define TK_ITHEORY 343
#define TK_NCOPIES 344
#define TK_NHITS 345
#define TK_MEAS_MESON_CORRS 346
#define TK_OPERATORS 347
#define TK_MEAS_NUCLEON_CORRS 348
#define TK_MEAS_PUTPOURRI 349
#define TK_COMPUTE_SUSC 350
#define TK_MAX_ORDER 351
#define TK_MEAS_RENDENS 352
#define TK_MEAS_ZUMBA 353
#define TK_MEAS_SPINPOL 354
#define TK_USE_FERM_CONF_FOR_GLUONS 355
#define TK_USE_ADJOINT_FLOW 356
#define TK_MEAS_QED_CORRS 357
#define TK_MEAS_MAGNETIZ 358
#define TK_MEAS_MIN_MAX_EIGENVAL 359
#define TK_MIN_MAX 360
#define TK_SPECTR_PROJ 361
#define TK_NEIGS 362
#define TK_EIG_PRECISION 363
#define TK_WSPACE_SIZE 364
#define TK_USE_SMOOTH 365
#define TK_MEAS_PLAQ_POL 366
#define TK_MEAS_PLAQ 367
#define TK_MEAS_ENERGY 368
#define TK_MEAS_POLY 369
#define TK_MEAS_TOP 370
#define TK_MEAS_LUPPOLI 371
#define TK_MEAS_WATUSSO 372
#define TK_MEAS_ALL_RECTS 373
#define TK_SPATIAL 374
#define TK_TEMPORAL 375
#define TK_DMIN 376
#define TK_DMAX 377
#define TK_TMIN 378
#define TK_TMAX 379
#define TK_EVOLUTION 380
#define TK_ID_SEA_THEORY 381
#define TK_NTRAJ_TOT 382
#define TK_SKIP_METRO 383
#define TK_TRAJ_LENGTH 384
#define TK_ACT_RESIDUE 385
#define TK_MD_RESIDUE 386
#define TK_NSUBSTEPS 387
#define TK_NPSEUDO_FERMS 388
#define TK_GAUGE_CONF 389
#define TK_STORE_PATH 390
#define TK_STORE_EACH 391
#define TK_STORE_RUNNING 392
#define TK_START_COND 393
#define TK_HOT 394
#define TK_COLD 395
#define TK_ANALYSIS 396
#define TK_CONF_LIST 397
#define TK_RUN 398
#define TK_WALLTIME 399
#define TK_SEED 400

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 37 "../../projects/generate_confs/parser.ypp" /* yacc.c:353  */

    double double_numb;
    int int_numb;
    std::string *text;
    nissa::gauge_action_name_t gauge_action_name;
    nissa::em_field_pars_t *em_field_pars;
    nissa::topotential_pars_t *topotential_pars;
    nissa::theory_pars_t *theory;
    nissa::quark_content_t *quark;
    nissa::smooth_pars_t *smooth_pars;
    nissa::stout_pars_t *stout_pars;
    nissa::cool_pars_t *cool_pars;
    nissa::ape_pars_t *ape_pars;
    nissa::hyp_pars_t *hyp_pars;
    nissa::Wflow_pars_t *Wflow_pars;
    
    nissa::meson_corr_meas_pars_t *meson_corr_meas;
    nissa::nucleon_corr_meas_pars_t *nucleon_corr_meas;
    nissa::magnetization_meas_pars_t *magnetization_meas;
    nissa::minmax_eigenvalues_meas_pars_t *minmax_eigenvalues_meas;
    nissa::quark_rendens_meas_pars_t *quark_rendens_meas;
    nissa::chir_zumba_meas_pars_t *chir_zumba_meas;
    nissa::spinpol_meas_pars_t *spinpol_meas;
    nissa::qed_corr_meas_pars_t *qed_corr_meas;
    nissa::fermionic_putpourri_meas_pars_t *fermionic_putpourri_meas;
    nissa::gauge_obs_meas_pars_t *plaq_pol_meas;
    nissa::top_meas_pars_t *top_meas;
    nissa::poly_corr_meas_pars_t *luppoli_meas;
    nissa::watusso_meas_pars_t *watusso_meas;
    nissa::all_rects_meas_pars_t *all_rects_meas;
    nissa::spectr_proj_meas_pars_t *spectr_proj_meas;
    
    nissa::rnd_t rnd_type;
    
    std::vector<int> *int_list;
    std::vector<double> *double_list;
    std::vector<std::string> *text_list;
    std::pair<int,int> *int_pair;
    std::vector<std::pair<int,int> > *int_pair_list;

#line 476 "generate_confs/parser.cpp" /* yacc.c:353  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif

/* Location type.  */
#if ! defined YYLTYPE && ! defined YYLTYPE_IS_DECLARED
typedef struct YYLTYPE YYLTYPE;
struct YYLTYPE
{
  int first_line;
  int first_column;
  int last_line;
  int last_column;
};
# define YYLTYPE_IS_DECLARED 1
# define YYLTYPE_IS_TRIVIAL 1
#endif



int parser_parse (nissa::driver_t *driver);

#endif /* !YY_PARSER_GENERATE_CONFS_PARSER_HPP_INCLUDED  */



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
typedef unsigned short yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short yytype_int16;
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
#  define YYSIZE_T unsigned
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

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && ! defined __ICC && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
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


#if 1

/* The parser invokes alloca or malloc; define the necessary symbols.  */

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
# define YYCOPY_NEEDED 1
#endif


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL \
             && defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
  YYLTYPE yyls_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE) + sizeof (YYLTYPE)) \
      + 2 * YYSTACK_GAP_MAXIMUM)

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
#define YYFINAL  56
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   645

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  153
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  74
/* YYNRULES -- Number of rules.  */
#define YYNRULES  322
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  539

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   400

#define YYTRANSLATE(YYX)                                                \
  ((unsigned) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     147,   148,     2,   152,   151,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   146,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   150,     2,   149,     2,     2,     2,     2,
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
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   126,   127,   128,   129,   130,   131,   132,   133,   134,
     135,   136,   137,   138,   139,   140,   141,   142,   143,   144,
     145
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   318,   318,   319,   320,   323,   324,   327,   329,   331,
     333,   334,   335,   336,   337,   338,   339,   340,   341,   343,
     344,   345,   346,   347,   348,   350,   351,   353,   355,   358,
     359,   360,   365,   366,   367,   369,   371,   373,   375,   377,
     379,   381,   383,   385,   387,   389,   391,   393,   395,   397,
     399,   400,   401,   402,   403,   404,   409,   410,   411,   412,
     413,   414,   415,   419,   420,   421,   426,   427,   428,   433,
     434,   435,   436,   437,   438,   439,   440,   441,   442,   443,
     444,   445,   450,   451,   452,   453,   454,   455,   456,   457,
     458,   463,   464,   465,   468,   469,   473,   474,   475,   478,
     482,   483,   484,   487,   491,   492,   493,   503,   507,   508,
     509,   510,   513,   514,   515,   519,   525,   526,   527,   528,
     529,   530,   531,   532,   533,   534,   535,   540,   541,   542,
     543,   544,   545,   546,   547,   548,   553,   554,   555,   556,
     557,   558,   559,   560,   561,   562,   567,   581,   582,   583,
     584,   585,   586,   587,   588,   589,   590,   595,   596,   597,
     598,   599,   600,   601,   602,   603,   604,   609,   610,   611,
     612,   613,   614,   615,   616,   617,   618,   623,   624,   625,
     626,   627,   628,   629,   630,   631,   632,   633,   634,   635,
     638,   639,   643,   644,   645,   646,   647,   648,   649,   650,
     651,   656,   657,   658,   659,   660,   661,   662,   663,   664,
     669,   670,   671,   672,   673,   674,   675,   676,   677,   678,
     679,   680,   685,   686,   687,   688,   689,   690,   691,   692,
     693,   694,   699,   700,   701,   711,   712,   713,   714,   719,
     720,   721,   722,   723,   724,   725,   726,   727,   728,   729,
     730,   735,   754,   755,   756,   757,   758,   759,   760,   761,
     762,   767,   768,   769,   770,   771,   772,   773,   778,   779,
     780,   781,   786,   787,   788,   789,   790,   791,   796,   797,
     798,   799,   800,   801,   802,   803,   804,   805,   810,   813,
     814,   819,   822,   823,   828,   831,   832,   837,   840,   841,
     846,   852,   853,   857,   858,   859,   860,   861,   862,   863,
     864,   865,   866,   870,   871,   872,   873,   874,   875,   876,
     877,   878,   879
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 1
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "TK_INVALID_STRING", "TK_INVALID_CHAR",
  "TK_DOUBLE", "TK_INT", "TK_QUANTITY", "TK_MINUS", "TK_PLUS", "TK_TIMES",
  "TK_DIV", "NEG", "TK_POW", "TK_GEOMETRY", "TK_LX", "TK_LY", "TK_LZ",
  "TK_L", "TK_T", "TK_TAG", "TK_THEORY", "TK_GAUGE_PARS", "TK_BETA",
  "TK_QUOTED_TEXT", "TK_EACH", "TK_MEAS_EACH_NSMOOTH", "TK_AFTER",
  "TK_NONE", "TK_YES", "TK_NO", "TK_RESIDUE", "TK_META", "TK_ORDINARY",
  "TK_TOPO_POTENTIAL", "TK_THETA", "TK_COEFF", "TK_WIDTH", "TK_BARR",
  "TK_FORCE_OUT", "TK_WELL_TEMPERING", "TK_BEND", "TK_PATH",
  "TK_META_PRINT_SCANNED_INPUT", "TK_META_PRINT_FULL_INPUT",
  "TK_META_PRINT", "TK_QUARK", "TK_DEGENERACY", "TK_DISCRETIZ",
  "TK_ROOT_STAG", "TK_OVERLAP", "TK_ROOT_TM_CLOV", "TK_KAPPA", "TK_CSW",
  "TK_MASS", "TK_MASS_OVERLAP", "TK_RE_POT_CH", "TK_IM_POT_CH",
  "TK_ELEC_CHARGE", "TK_SMOOTH_METHOD", "TK_SPACE_OR_TIME", "TK_SPACE",
  "TK_TIME", "TK_SPACETIME", "TK_MEAS_CORR", "TK_CORR_PATH",
  "TK_NOISE_TYPE", "TK_RND_T", "TK_STOUT", "TK_NLEVELS", "TK_RHO",
  "TK_COOLING", "TK_NSTEPS", "TK_APE", "TK_ALPHA", "TK_HYP", "TK_ALPHAS",
  "TK_WFLOW", "TK_NFLOWS", "TK_NRECU", "TK_FLOW_STEP", "TK_GAUGE_ACTION",
  "TK_WILSON", "TK_TLSYM", "TK_IWASAKI", "TK_BKGRD_EM_FIELD", "TK_B_COMP",
  "TK_E_COMP", "TK_ITHEORY", "TK_NCOPIES", "TK_NHITS",
  "TK_MEAS_MESON_CORRS", "TK_OPERATORS", "TK_MEAS_NUCLEON_CORRS",
  "TK_MEAS_PUTPOURRI", "TK_COMPUTE_SUSC", "TK_MAX_ORDER",
  "TK_MEAS_RENDENS", "TK_MEAS_ZUMBA", "TK_MEAS_SPINPOL",
  "TK_USE_FERM_CONF_FOR_GLUONS", "TK_USE_ADJOINT_FLOW",
  "TK_MEAS_QED_CORRS", "TK_MEAS_MAGNETIZ", "TK_MEAS_MIN_MAX_EIGENVAL",
  "TK_MIN_MAX", "TK_SPECTR_PROJ", "TK_NEIGS", "TK_EIG_PRECISION",
  "TK_WSPACE_SIZE", "TK_USE_SMOOTH", "TK_MEAS_PLAQ_POL", "TK_MEAS_PLAQ",
  "TK_MEAS_ENERGY", "TK_MEAS_POLY", "TK_MEAS_TOP", "TK_MEAS_LUPPOLI",
  "TK_MEAS_WATUSSO", "TK_MEAS_ALL_RECTS", "TK_SPATIAL", "TK_TEMPORAL",
  "TK_DMIN", "TK_DMAX", "TK_TMIN", "TK_TMAX", "TK_EVOLUTION",
  "TK_ID_SEA_THEORY", "TK_NTRAJ_TOT", "TK_SKIP_METRO", "TK_TRAJ_LENGTH",
  "TK_ACT_RESIDUE", "TK_MD_RESIDUE", "TK_NSUBSTEPS", "TK_NPSEUDO_FERMS",
  "TK_GAUGE_CONF", "TK_STORE_PATH", "TK_STORE_EACH", "TK_STORE_RUNNING",
  "TK_START_COND", "TK_HOT", "TK_COLD", "TK_ANALYSIS", "TK_CONF_LIST",
  "TK_RUN", "TK_WALLTIME", "TK_SEED", "'='", "'('", "')'", "'}'", "'{'",
  "','", "'+'", "$accept", "commands", "command", "global_specify",
  "meta_command", "gauge_action", "each", "meas_each_nsmooth", "after",
  "residue", "itheory", "ncopies", "nhits", "noise_type", "compute_susc",
  "neigs", "eig_precision", "wspace_size", "min_max", "max_order", "path",
  "geometry", "theory", "em_field_pars", "run_pars", "topo_potential_pars",
  "smooth_pars", "stout_pars", "nlevels", "rho", "cool_pars", "nsteps",
  "ape_pars", "alpha", "hyp_pars", "alphas", "Wflow_pars", "nflows",
  "nrecu", "flow_step", "quark", "nucleon_corr_meas", "meson_corr_meas",
  "operators", "fermionic_putpourri_meas", "quark_rendens_meas",
  "chir_zumba_meas", "spinpol_meas", "use_ferm_conf_for_gluons",
  "use_adjoint_flow", "qed_corr_meas", "magnetization_meas",
  "minmax_eigenvalues_meas", "global_evolve_pars", "global_conf_pars",
  "spectr_proj_meas", "global_analysis_pars", "plaq_pol_meas", "top_meas",
  "luppoli_meas", "watusso_meas", "all_rects_meas", "int_list",
  "internal_int_list", "double_list", "internal_double_list", "text_list",
  "internal_text_list", "int_pair_list", "internal_int_pair_list",
  "int_pair", "text", "double_numb", "int_numb", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   379,   380,   381,   382,   383,   384,
     385,   386,   387,   388,   389,   390,   391,   392,   393,   394,
     395,   396,   397,   398,   399,   400,    61,    40,    41,   125,
     123,    44,    43
};
# endif

#define YYPACT_NINF -237

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-237)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     340,  -237,  -124,  -237,  -105,   -83,   -67,  -237,  -237,  -237,
    -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,
    -237,  -237,  -237,  -237,   -48,  -237,   253,  -237,  -237,  -237,
     199,    87,   -73,   473,   437,   277,   245,   245,   442,   473,
     473,   382,    20,    -7,   387,  -237,    99,    73,    92,    40,
      76,    75,    75,    75,    75,   -44,  -237,  -237,   -23,   -11,
       8,    21,    28,    30,    47,    55,  -237,    64,  -237,  -237,
     121,   285,   155,   464,    81,   117,   120,   122,   123,   136,
     144,   146,   148,   153,  -237,  -237,  -237,  -237,  -237,  -237,
    -237,  -237,   154,  -237,  -237,  -237,  -237,  -237,  -237,  -237,
    -237,  -237,   157,  -237,  -237,  -237,  -237,  -237,  -237,  -237,
    -237,  -237,   159,  -237,  -237,  -237,  -237,  -237,  -237,  -237,
    -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,
    -237,   160,   161,   163,  -237,  -237,  -237,  -237,  -237,  -237,
    -237,  -237,   -17,  -237,  -237,  -237,  -237,  -237,  -237,  -237,
    -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,
    -237,  -237,   171,   181,   202,  -237,  -237,  -237,  -237,  -237,
    -237,  -237,  -237,  -237,  -237,  -237,   203,   212,   216,   227,
     228,   229,   230,   231,   233,   234,   236,   240,   242,  -237,
     243,  -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,
    -237,  -237,   244,   246,   247,   249,  -237,  -237,  -237,   -17,
     252,   254,  -237,  -237,  -237,   -17,  -237,  -237,  -237,   105,
     105,  -237,  -237,  -237,   105,   105,   256,   257,   258,   259,
    -237,  -237,  -237,  -237,    70,    25,    33,    34,   133,    24,
      24,    24,    24,    24,    15,   137,    75,    61,   262,   264,
     265,   270,   271,   273,   274,   275,   276,  -237,  -237,   155,
     279,   280,  -237,  -237,   281,   282,   284,   286,   289,   290,
     294,   295,   299,    24,    24,    24,    24,    15,    75,   278,
      24,    24,    24,   151,    24,    24,   218,    24,    24,   301,
     303,  -237,    24,    24,    24,    24,    24,    24,    24,    15,
      15,    15,    24,   213,    75,    24,    24,    90,    15,    24,
      24,    24,    24,    24,    75,   -17,   -17,   -17,   -17,    24,
      24,    24,    24,    75,  -237,  -237,  -237,    75,  -237,  -115,
    -237,    24,    24,    24,   329,   329,   329,   329,   329,  -237,
      15,    15,    15,   305,   329,  -237,  -237,  -237,    70,  -237,
    -237,  -237,    15,    15,    15,    15,    15,    15,    15,    15,
      15,    24,    15,    24,    58,    15,    15,    15,    15,    15,
      15,    15,   329,   329,   329,   329,   305,    70,  -237,   329,
     329,   329,   307,  -237,   -96,   329,   329,  -237,  -237,  -237,
    -237,   155,   -31,    -6,   -22,   112,   329,   329,    15,   217,
     329,   329,   329,   329,   329,   329,   329,   305,   305,   305,
     329,    24,  -237,   -29,    70,   329,   329,  -237,  -237,   305,
     329,   329,   329,   329,   329,    70,   329,   329,   329,   329,
      70,    70,  -237,    75,     5,     5,    18,  -237,    24,    24,
      24,    24,    24,   378,     5,   378,     5,   103,    18,    15,
      15,    15,    15,    15,   305,   305,   305,   305,   305,   305,
     305,   305,   305,   329,   305,   329,  -237,  -237,  -237,   305,
     305,   305,   305,   305,   305,   305,    24,  -237,  -237,   307,
     304,  -237,  -237,   306,  -237,  -237,   313,  -237,  -237,   314,
     315,   317,  -237,  -237,  -237,   305,  -237,  -237,  -237,   329,
    -237,    24,    70,  -237,   129,   129,     5,     5,     5,  -237,
      68,    68,   378,   378,   378,     6,  -237,    24,    15,   316,
      24,    24,    15,   329,    24,   329,   305,    15,  -237,    12,
     329,   329,   305,    80,   305,  -237,    15,  -237,   305
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint16 yydefact[] =
{
       2,    50,     0,    56,     0,     0,     0,   136,   127,   147,
     157,   167,   177,   192,   201,   210,   239,   252,   261,   268,
     272,   278,   222,   232,     0,    66,     0,     4,     5,     6,
       7,     9,    28,    11,    10,    12,    13,    14,    15,    16,
      17,    18,    25,    26,    21,    27,    19,    20,    22,    23,
      24,     0,     0,     0,     0,     0,     1,     3,     0,     0,
       0,     0,     0,     0,     0,     0,    91,     0,    63,    62,
      60,    58,    57,    59,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   128,   129,   131,   132,   133,   135,
     134,   130,     0,   137,   138,   140,   141,   142,   144,   143,
     139,   145,     0,   148,   149,   151,   153,   154,   156,   155,
     152,   150,     0,   158,   159,   161,   162,   163,   165,   164,
     166,   160,   168,   169,   171,   172,   173,   175,   174,   176,
     170,     0,     0,     0,   178,   179,   181,   182,   183,   185,
     184,   180,   189,   186,   187,   188,   193,   194,   196,   197,
     198,   200,   199,   195,   202,   203,   205,   206,   207,   209,
     208,   204,     0,     0,     0,   211,   212,   214,   215,   216,
     218,   217,   219,   220,   221,   213,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   233,
       0,   240,   241,   243,   244,   245,   247,   246,   248,   249,
     250,   242,     0,     0,     0,     0,   253,   254,   255,   260,
       0,     0,   262,   263,   264,   265,   269,   270,   271,     0,
       0,   273,   274,   277,     0,     0,     0,     0,     0,     0,
     279,   280,   287,   301,     8,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    74,    73,    72,
       0,     0,    92,    93,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    87,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   275,   276,   285,   286,     0,
       0,     0,     0,     0,    29,    30,    31,     0,   251,     0,
     313,     0,     0,     0,    53,    54,    55,    52,    51,   303,
       0,     0,     0,    61,   304,    71,    69,    70,   115,    32,
      33,    34,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    67,    68,    35,    37,    38,    49,    42,    39,
      40,    41,     0,   146,     0,    43,    48,    96,   100,   104,
     108,    82,    83,    85,    86,    84,   190,   191,     0,     0,
      47,    44,    46,   229,   224,   223,   225,   226,   227,   228,
     230,     0,   231,     0,   234,   235,   236,   237,   238,    45,
     256,   257,   258,   259,   266,   267,   281,   282,   283,   284,
     302,   295,   294,     0,   315,   314,     0,   318,     0,     0,
       0,     0,     0,   306,   304,   305,   304,     0,     0,     0,
       0,     0,     0,     0,    65,    64,    75,    76,    77,    78,
      79,    80,    81,    94,    95,   116,   117,   119,   118,   122,
     123,   120,   121,   124,   125,   126,     0,   298,   297,     0,
       0,    98,    97,     0,   101,   102,     0,   105,   106,     0,
       0,     0,   109,   110,   111,    36,    88,    89,    90,   289,
     288,     0,   296,   322,   317,   316,   319,   320,   321,   312,
     308,   307,   309,   310,   311,     0,   299,     0,     0,     0,
       0,     0,     0,   290,     0,    99,   103,     0,   107,     0,
     112,   113,   114,     0,   292,   291,     0,   300,   293
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -237,  -237,   373,  -237,  -237,    31,   510,  -237,   530,   220,
     568,   579,   590,   601,  -237,   434,  -237,   436,  -237,   445,
     548,  -237,  -237,  -237,  -237,  -237,   -36,   -63,  -151,  -237,
    -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,
    -237,  -237,  -237,   447,  -237,  -237,  -237,  -237,  -237,  -237,
    -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,
    -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,  -237,
       7,   -52,  -121,  -236
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    26,    27,    28,    29,    69,    84,   291,    85,    86,
      87,    88,    89,    90,   110,   172,   199,   173,   174,   120,
      91,    30,    31,    70,    32,    71,   142,    72,   262,   263,
     392,   482,   393,   485,   394,   488,   395,   492,   493,   494,
      73,    33,    34,   101,    35,    36,    37,    38,   144,   145,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,   412,   413,   528,   529,   328,   329,   383,   384,
     477,   234,   343,   344
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_uint16 yytable[] =
{
     235,   236,   237,   334,   335,   336,   337,   338,   259,   289,
     209,   215,   437,   437,   438,   439,   440,   441,   442,   442,
     339,   330,    51,   340,   341,   437,   438,   439,   440,   441,
     330,   442,   331,   332,   432,    79,   433,   372,   373,   374,
     375,   480,    52,   290,   379,   380,   381,   260,   385,   386,
      67,   396,   397,   478,   486,   479,   400,   401,   402,   403,
     404,   405,   406,   260,    53,    76,   410,    77,   483,   415,
     416,    74,    75,   420,   421,   422,   423,   424,   451,   452,
      54,   453,    79,   426,   427,   428,   429,   437,   438,   439,
     440,   441,   176,   442,    55,   434,   435,   436,    76,   233,
      77,    76,   238,    77,   444,   446,   448,   466,   467,   468,
      63,   449,   450,   451,   452,    79,   453,    76,    79,    77,
     500,    64,   501,   239,    76,   463,    77,   465,   185,   186,
     187,   188,   131,    65,    79,   240,   437,   210,   211,   440,
     441,    79,   442,   349,   350,   351,   177,   178,   179,   180,
     181,   182,   183,   184,   241,    66,   376,   524,   131,   219,
     220,   535,   342,   536,   131,   345,   503,   242,    67,   346,
     347,   333,    68,   324,   243,   499,   244,   323,   407,   408,
     409,   325,   326,   315,   316,   323,   323,   419,   317,   318,
     489,   490,   491,   245,   348,   224,   225,   226,   227,   228,
     229,   246,   504,   505,   506,   507,   508,   248,   249,   202,
     247,   203,   204,   205,    58,    59,    60,    61,    62,   443,
     445,   447,   323,   391,   260,   261,   377,   273,   537,   417,
     418,   454,   455,   456,   457,   458,   459,   460,   461,   462,
     515,   464,   484,   487,   469,   470,   471,   472,   473,   474,
     475,   509,   414,    56,    95,   105,   115,   124,   136,   148,
     156,   167,   425,   274,   193,   523,   275,     1,   276,   277,
      76,   430,    77,     2,     3,   431,    78,   495,   496,   497,
     498,   525,   278,   327,   530,   531,    66,    79,   533,   387,
     279,   388,   280,   389,   281,   390,     4,     5,     6,   282,
     283,   382,    76,   284,    77,   285,   286,   287,    78,   288,
      76,    80,    77,   449,   450,   451,   452,   292,   453,    79,
     250,   251,   252,   253,   254,   255,   256,   293,   510,   511,
     512,   513,   514,    81,    82,    83,   437,   438,   439,   440,
     441,   112,   442,    80,     7,   378,     8,     9,   294,   295,
      10,    11,    12,    66,     1,    13,    14,    15,   296,    16,
       2,     3,   297,   411,    17,    81,    82,    83,    18,    19,
      20,    21,   102,   298,   299,   300,   301,   302,    22,   303,
     304,   502,   305,     4,     5,     6,   306,    23,   307,   308,
     309,   453,   310,   311,    24,   312,    25,   526,   313,    57,
     314,   532,   319,   320,   321,   322,   534,    76,   352,    77,
     353,   354,    76,    78,    77,   538,   355,   356,    78,   357,
     358,   359,   360,   481,    79,   361,   362,   363,   364,    79,
     365,     7,   366,     8,     9,   367,   368,    10,    11,    12,
     369,   370,    13,    14,    15,   371,    16,   398,    80,   399,
     517,    17,   518,    80,   476,    18,    19,    20,    21,   519,
     520,   521,    76,   522,    77,    22,   527,    76,    78,    77,
      81,    82,    83,    78,    23,    81,    82,    83,   198,    79,
     200,    24,   129,    25,    79,   143,   516,   162,     0,   163,
       0,   164,     0,     0,   163,   190,   164,     0,    76,     0,
      77,   131,     0,    80,    78,     0,     0,     0,    80,     0,
       0,   264,   265,     0,     0,    79,   266,   267,   268,   269,
     270,   271,   272,     0,     0,    81,    82,    83,     0,    92,
      81,    82,    83,     0,    92,     0,     0,     0,     0,    80,
       0,     0,   132,   133,    93,   103,   113,   122,   134,   146,
     154,   165,     0,     0,   191,     0,   206,   212,   216,   221,
     230,    81,    82,    83,    94,   104,   114,   123,   135,   147,
     155,   166,     0,     0,   192,     0,   207,   213,   217,   222,
     231,   257,   100,   111,   121,   130,   141,   153,   161,   175,
       0,   189,   201,     0,   208,   214,   218,   223,   232,     0,
       0,   258,    96,   106,   116,   125,   137,   149,   157,   168,
       0,     0,   194,    97,   107,   117,   126,   138,   150,   158,
     169,     0,     0,   195,    98,   108,   118,   127,   139,   151,
     159,   170,     0,     0,   196,    99,   109,   119,   128,   140,
     152,   160,   171,     0,     0,   197
};

static const yytype_int16 yycheck[] =
{
      52,    53,    54,   239,   240,   241,   242,   243,    71,    26,
      46,    47,     7,     7,     8,     9,    10,    11,    13,    13,
       5,     6,   146,     8,     9,     7,     8,     9,    10,    11,
       6,    13,     8,     9,   149,    42,   151,   273,   274,   275,
     276,    72,   147,    60,   280,   281,   282,    69,   284,   285,
      81,   287,   288,   149,    76,   151,   292,   293,   294,   295,
     296,   297,   298,    69,   147,    25,   302,    27,    74,   305,
     306,   144,   145,   309,   310,   311,   312,   313,    10,    11,
     147,    13,    42,   319,   320,   321,   322,     7,     8,     9,
      10,    11,    72,    13,   142,   331,   332,   333,    25,    24,
      27,    25,   146,    27,   340,   341,   342,    49,    50,    51,
      23,     8,     9,    10,    11,    42,    13,    25,    42,    27,
     149,    34,   151,   146,    25,   361,    27,   363,   135,   136,
     137,   138,    59,    46,    42,   146,     7,    64,    65,    10,
      11,    42,    13,    82,    83,    84,   126,   127,   128,   129,
     130,   131,   132,   133,   146,    68,   277,   151,    59,   119,
     120,   149,   147,   151,    59,    28,   148,   146,    81,    32,
      33,   147,    85,   148,   146,   411,   146,   152,   299,   300,
     301,   148,   148,   219,   220,   152,   152,   308,   224,   225,
      78,    79,    80,   146,   246,   119,   120,   121,   122,   123,
     124,   146,   438,   439,   440,   441,   442,    86,    87,   110,
     146,   112,   113,   114,    15,    16,    17,    18,    19,   340,
     341,   342,   152,   286,    69,    70,   278,   146,   148,   139,
     140,   352,   353,   354,   355,   356,   357,   358,   359,   360,
     476,   362,   393,   394,   365,   366,   367,   368,   369,   370,
     371,   148,   304,     0,    34,    35,    36,    37,    38,    39,
      40,    41,   314,   146,    44,   501,   146,    14,   146,   146,
      25,   323,    27,    20,    21,   327,    31,   398,    61,    62,
      63,   517,   146,   150,   520,   521,    68,    42,   524,    71,
     146,    73,   146,    75,   146,    77,    43,    44,    45,   146,
     146,   150,    25,   146,    27,   146,   146,   146,    31,   146,
      25,    66,    27,     8,     9,    10,    11,   146,    13,    42,
      35,    36,    37,    38,    39,    40,    41,   146,   449,   450,
     451,   452,   453,    88,    89,    90,     7,     8,     9,    10,
      11,    96,    13,    66,    91,    67,    93,    94,   146,   146,
      97,    98,    99,    68,    14,   102,   103,   104,   146,   106,
      20,    21,   146,   150,   111,    88,    89,    90,   115,   116,
     117,   118,    95,   146,   146,   146,   146,   146,   125,   146,
     146,   433,   146,    43,    44,    45,   146,   134,   146,   146,
     146,    13,   146,   146,   141,   146,   143,   518,   146,    26,
     146,   522,   146,   146,   146,   146,   527,    25,   146,    27,
     146,   146,    25,    31,    27,   536,   146,   146,    31,   146,
     146,   146,   146,   392,    42,   146,   146,   146,   146,    42,
     146,    91,   146,    93,    94,   146,   146,    97,    98,    99,
     146,   146,   102,   103,   104,   146,   106,   146,    66,   146,
     146,   111,   146,    66,   147,   115,   116,   117,   118,   146,
     146,   146,    25,   146,    27,   125,   150,    25,    31,    27,
      88,    89,    90,    31,   134,    88,    89,    90,    44,    42,
      44,   141,    37,   143,    42,    38,   479,   105,    -1,   107,
      -1,   109,    -1,    -1,   107,   108,   109,    -1,    25,    -1,
      27,    59,    -1,    66,    31,    -1,    -1,    -1,    66,    -1,
      -1,    47,    48,    -1,    -1,    42,    52,    53,    54,    55,
      56,    57,    58,    -1,    -1,    88,    89,    90,    -1,    92,
      88,    89,    90,    -1,    92,    -1,    -1,    -1,    -1,    66,
      -1,    -1,   100,   101,    34,    35,    36,    37,    38,    39,
      40,    41,    -1,    -1,    44,    -1,    46,    47,    48,    49,
      50,    88,    89,    90,    34,    35,    36,    37,    38,    39,
      40,    41,    -1,    -1,    44,    -1,    46,    47,    48,    49,
      50,    71,    34,    35,    36,    37,    38,    39,    40,    41,
      -1,    43,    44,    -1,    46,    47,    48,    49,    50,    -1,
      -1,    71,    34,    35,    36,    37,    38,    39,    40,    41,
      -1,    -1,    44,    34,    35,    36,    37,    38,    39,    40,
      41,    -1,    -1,    44,    34,    35,    36,    37,    38,    39,
      40,    41,    -1,    -1,    44,    34,    35,    36,    37,    38,
      39,    40,    41,    -1,    -1,    44
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    14,    20,    21,    43,    44,    45,    91,    93,    94,
      97,    98,    99,   102,   103,   104,   106,   111,   115,   116,
     117,   118,   125,   134,   141,   143,   154,   155,   156,   157,
     174,   175,   177,   194,   195,   197,   198,   199,   200,   203,
     204,   205,   206,   207,   208,   209,   210,   211,   212,   213,
     214,   146,   147,   147,   147,   142,     0,   155,    15,    16,
      17,    18,    19,    23,    34,    46,    68,    81,    85,   158,
     176,   178,   180,   193,   144,   145,    25,    27,    31,    42,
      66,    88,    89,    90,   159,   161,   162,   163,   164,   165,
     166,   173,    92,   159,   161,   162,   163,   164,   165,   166,
     173,   196,    95,   159,   161,   162,   163,   164,   165,   166,
     167,   173,    96,   159,   161,   162,   163,   164,   165,   166,
     172,   173,   159,   161,   162,   163,   164,   165,   166,   172,
     173,    59,   100,   101,   159,   161,   162,   163,   164,   165,
     166,   173,   179,   196,   201,   202,   159,   161,   162,   163,
     164,   165,   166,   173,   159,   161,   162,   163,   164,   165,
     166,   173,   105,   107,   109,   159,   161,   162,   163,   164,
     165,   166,   168,   170,   171,   173,    72,   126,   127,   128,
     129,   130,   131,   132,   133,   135,   136,   137,   138,   173,
     108,   159,   161,   162,   163,   164,   165,   166,   168,   169,
     170,   173,   110,   112,   113,   114,   159,   161,   173,   179,
      64,    65,   159,   161,   173,   179,   159,   161,   173,   119,
     120,   159,   161,   173,   119,   120,   121,   122,   123,   124,
     159,   161,   173,    24,   224,   224,   224,   224,   146,   146,
     146,   146,   146,   146,   146,   146,   146,   146,    86,    87,
      35,    36,    37,    38,    39,    40,    41,   159,   161,   180,
      69,    70,   181,   182,    47,    48,    52,    53,    54,    55,
      56,    57,    58,   146,   146,   146,   146,   146,   146,   146,
     146,   146,   146,   146,   146,   146,   146,   146,   146,    26,
      60,   160,   146,   146,   146,   146,   146,   146,   146,   146,
     146,   146,   146,   146,   146,   146,   146,   146,   146,   146,
     146,   146,   146,   146,   146,   179,   179,   179,   179,   146,
     146,   146,   146,   152,   148,   148,   148,   150,   219,   220,
       6,     8,     9,   147,   226,   226,   226,   226,   226,     5,
       8,     9,   147,   225,   226,    28,    32,    33,   224,    82,
      83,    84,   146,   146,   146,   146,   146,   146,   146,   146,
     146,   146,   146,   146,   146,   146,   146,   146,   146,   146,
     146,   146,   226,   226,   226,   226,   225,   224,    67,   226,
     226,   226,   150,   221,   222,   226,   226,    71,    73,    75,
      77,   180,   183,   185,   187,   189,   226,   226,   146,   146,
     226,   226,   226,   226,   226,   226,   226,   225,   225,   225,
     226,   150,   215,   216,   224,   226,   226,   139,   140,   225,
     226,   226,   226,   226,   226,   224,   226,   226,   226,   226,
     224,   224,   149,   151,   226,   226,   226,     7,     8,     9,
      10,    11,    13,   225,   226,   225,   226,   225,   226,     8,
       9,    10,    11,    13,   225,   225,   225,   225,   225,   225,
     225,   225,   225,   226,   225,   226,    49,    50,    51,   225,
     225,   225,   225,   225,   225,   225,   147,   223,   149,   151,
      72,   158,   184,    74,   181,   186,    76,   181,   188,    78,
      79,    80,   190,   191,   192,   225,    61,    62,    63,   226,
     149,   151,   224,   148,   226,   226,   226,   226,   226,   148,
     225,   225,   225,   225,   225,   226,   223,   146,   146,   146,
     146,   146,   146,   226,   151,   226,   225,   150,   217,   218,
     226,   226,   225,   226,   225,   149,   151,   148,   225
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,   153,   154,   154,   154,   155,   155,   156,   156,   156,
     156,   156,   156,   156,   156,   156,   156,   156,   156,   156,
     156,   156,   156,   156,   156,   156,   156,   156,   156,   157,
     157,   157,   158,   158,   158,   159,   160,   161,   162,   163,
     164,   165,   166,   167,   168,   169,   170,   171,   172,   173,
     174,   174,   174,   174,   174,   174,   175,   175,   175,   175,
     175,   175,   175,   176,   176,   176,   177,   177,   177,   178,
     178,   178,   178,   178,   178,   178,   178,   178,   178,   178,
     178,   178,   179,   179,   179,   179,   179,   179,   179,   179,
     179,   180,   180,   180,   181,   182,   183,   183,   183,   184,
     185,   185,   185,   186,   187,   187,   187,   188,   189,   189,
     189,   189,   190,   191,   192,   193,   193,   193,   193,   193,
     193,   193,   193,   193,   193,   193,   193,   194,   194,   194,
     194,   194,   194,   194,   194,   194,   195,   195,   195,   195,
     195,   195,   195,   195,   195,   195,   196,   197,   197,   197,
     197,   197,   197,   197,   197,   197,   197,   198,   198,   198,
     198,   198,   198,   198,   198,   198,   198,   199,   199,   199,
     199,   199,   199,   199,   199,   199,   199,   200,   200,   200,
     200,   200,   200,   200,   200,   200,   200,   200,   200,   200,
     201,   202,   203,   203,   203,   203,   203,   203,   203,   203,
     203,   204,   204,   204,   204,   204,   204,   204,   204,   204,
     205,   205,   205,   205,   205,   205,   205,   205,   205,   205,
     205,   205,   206,   206,   206,   206,   206,   206,   206,   206,
     206,   206,   207,   207,   207,   207,   207,   207,   207,   208,
     208,   208,   208,   208,   208,   208,   208,   208,   208,   208,
     208,   209,   210,   210,   210,   210,   210,   210,   210,   210,
     210,   211,   211,   211,   211,   211,   211,   211,   212,   212,
     212,   212,   213,   213,   213,   213,   213,   213,   214,   214,
     214,   214,   214,   214,   214,   214,   214,   214,   215,   216,
     216,   217,   218,   218,   219,   220,   220,   221,   222,   222,
     223,   224,   224,   225,   225,   225,   225,   225,   225,   225,
     225,   225,   225,   226,   226,   226,   226,   226,   226,   226,
     226,   226,   226
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     1,     1,     1,     1,     3,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     4,
       4,     4,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       1,     4,     4,     4,     4,     4,     1,     2,     2,     2,
       2,     4,     2,     1,     4,     4,     1,     4,     4,     3,
       3,     3,     2,     2,     2,     4,     4,     4,     4,     4,
       4,     4,     3,     3,     3,     3,     3,     2,     4,     4,
       4,     1,     2,     2,     3,     3,     1,     2,     2,     3,
       1,     2,     2,     3,     1,     2,     2,     3,     1,     2,
       2,     2,     3,     3,     3,     3,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     1,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     3,     1,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     1,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     1,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     1,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       3,     3,     1,     2,     2,     2,     2,     2,     2,     2,
       2,     1,     2,     2,     2,     2,     2,     2,     2,     2,
       1,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     1,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     1,     2,     4,     4,     4,     4,     4,     1,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     4,     1,     2,     2,     2,     4,     4,     4,     4,
       2,     1,     2,     2,     2,     2,     4,     4,     1,     2,
       2,     2,     1,     2,     2,     3,     3,     2,     1,     2,
       2,     4,     4,     4,     4,     3,     3,     2,     2,     2,
       3,     2,     2,     3,     2,     2,     3,     2,     2,     3,
       5,     1,     3,     1,     1,     2,     2,     3,     3,     3,
       3,     3,     3,     1,     2,     2,     3,     3,     2,     3,
       3,     3,     3
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
      YY_LAC_DISCARD ("YYBACKUP");                              \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (&yylloc, driver, YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)                                \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;        \
          (Current).first_column = YYRHSLOC (Rhs, 1).first_column;      \
          (Current).last_line    = YYRHSLOC (Rhs, N).last_line;         \
          (Current).last_column  = YYRHSLOC (Rhs, N).last_column;       \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).first_line   = (Current).last_line   =              \
            YYRHSLOC (Rhs, 0).last_line;                                \
          (Current).first_column = (Current).last_column =              \
            YYRHSLOC (Rhs, 0).last_column;                              \
        }                                                               \
    while (0)
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K])


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


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL

/* Print *YYLOCP on YYO.  Private, do not rely on its existence. */

YY_ATTRIBUTE_UNUSED
static int
yy_location_print_ (FILE *yyo, YYLTYPE const * const yylocp)
{
  int res = 0;
  int end_col = 0 != yylocp->last_column ? yylocp->last_column - 1 : 0;
  if (0 <= yylocp->first_line)
    {
      res += YYFPRINTF (yyo, "%d", yylocp->first_line);
      if (0 <= yylocp->first_column)
        res += YYFPRINTF (yyo, ".%d", yylocp->first_column);
    }
  if (0 <= yylocp->last_line)
    {
      if (yylocp->first_line < yylocp->last_line)
        {
          res += YYFPRINTF (yyo, "-%d", yylocp->last_line);
          if (0 <= end_col)
            res += YYFPRINTF (yyo, ".%d", end_col);
        }
      else if (0 <= end_col && yylocp->first_column < end_col)
        res += YYFPRINTF (yyo, "-%d", end_col);
    }
  return res;
 }

#  define YY_LOCATION_PRINT(File, Loc)          \
  yy_location_print_ (File, &(Loc))

# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value, Location, driver); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp, nissa::driver_t *driver)
{
  FILE *yyoutput = yyo;
  YYUSE (yyoutput);
  YYUSE (yylocationp);
  YYUSE (driver);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyo, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp, nissa::driver_t *driver)
{
  YYFPRINTF (yyo, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  YY_LOCATION_PRINT (yyo, *yylocationp);
  YYFPRINTF (yyo, ": ");
  yy_symbol_value_print (yyo, yytype, yyvaluep, yylocationp, driver);
  YYFPRINTF (yyo, ")");
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
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, YYLTYPE *yylsp, int yyrule, nissa::driver_t *driver)
{
  unsigned long yylno = yyrline[yyrule];
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
                       , &(yylsp[(yyi + 1) - (yynrhs)])                       , driver);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, yylsp, Rule, driver); \
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

/* Given a state stack such that *YYBOTTOM is its bottom, such that
   *YYTOP is either its top or is YYTOP_EMPTY to indicate an empty
   stack, and such that *YYCAPACITY is the maximum number of elements it
   can hold without a reallocation, make sure there is enough room to
   store YYADD more elements.  If not, allocate a new stack using
   YYSTACK_ALLOC, copy the existing elements, and adjust *YYBOTTOM,
   *YYTOP, and *YYCAPACITY to reflect the new capacity and memory
   location.  If *YYBOTTOM != YYBOTTOM_NO_FREE, then free the old stack
   using YYSTACK_FREE.  Return 0 if successful or if no reallocation is
   required.  Return 1 if memory is exhausted.  */
static int
yy_lac_stack_realloc (YYSIZE_T *yycapacity, YYSIZE_T yyadd,
#if YYDEBUG
                      char const *yydebug_prefix,
                      char const *yydebug_suffix,
#endif
                      yytype_int16 **yybottom,
                      yytype_int16 *yybottom_no_free,
                      yytype_int16 **yytop, yytype_int16 *yytop_empty)
{
  YYSIZE_T yysize_old =
    (YYSIZE_T) (*yytop == yytop_empty ? 0 : *yytop - *yybottom + 1);
  YYSIZE_T yysize_new = yysize_old + yyadd;
  if (*yycapacity < yysize_new)
    {
      YYSIZE_T yyalloc = 2 * yysize_new;
      yytype_int16 *yybottom_new;
      /* Use YYMAXDEPTH for maximum stack size given that the stack
         should never need to grow larger than the main state stack
         needs to grow without LAC.  */
      if (YYMAXDEPTH < yysize_new)
        {
          YYDPRINTF ((stderr, "%smax size exceeded%s", yydebug_prefix,
                      yydebug_suffix));
          return 1;
        }
      if (YYMAXDEPTH < yyalloc)
        yyalloc = YYMAXDEPTH;
      yybottom_new =
        (yytype_int16*) YYSTACK_ALLOC (yyalloc * sizeof *yybottom_new);
      if (!yybottom_new)
        {
          YYDPRINTF ((stderr, "%srealloc failed%s", yydebug_prefix,
                      yydebug_suffix));
          return 1;
        }
      if (*yytop != yytop_empty)
        {
          YYCOPY (yybottom_new, *yybottom, yysize_old);
          *yytop = yybottom_new + (yysize_old - 1);
        }
      if (*yybottom != yybottom_no_free)
        YYSTACK_FREE (*yybottom);
      *yybottom = yybottom_new;
      *yycapacity = yyalloc;
    }
  return 0;
}

/* Establish the initial context for the current lookahead if no initial
   context is currently established.

   We define a context as a snapshot of the parser stacks.  We define
   the initial context for a lookahead as the context in which the
   parser initially examines that lookahead in order to select a
   syntactic action.  Thus, if the lookahead eventually proves
   syntactically unacceptable (possibly in a later context reached via a
   series of reductions), the initial context can be used to determine
   the exact set of tokens that would be syntactically acceptable in the
   lookahead's place.  Moreover, it is the context after which any
   further semantic actions would be erroneous because they would be
   determined by a syntactically unacceptable token.

   YY_LAC_ESTABLISH should be invoked when a reduction is about to be
   performed in an inconsistent state (which, for the purposes of LAC,
   includes consistent states that don't know they're consistent because
   their default reductions have been disabled).  Iff there is a
   lookahead token, it should also be invoked before reporting a syntax
   error.  This latter case is for the sake of the debugging output.

   For parse.lac=full, the implementation of YY_LAC_ESTABLISH is as
   follows.  If no initial context is currently established for the
   current lookahead, then check if that lookahead can eventually be
   shifted if syntactic actions continue from the current context.
   Report a syntax error if it cannot.  */
#define YY_LAC_ESTABLISH                                         \
do {                                                             \
  if (!yy_lac_established)                                       \
    {                                                            \
      YYDPRINTF ((stderr,                                        \
                  "LAC: initial context established for %s\n",   \
                  yytname[yytoken]));                            \
      yy_lac_established = 1;                                    \
      {                                                          \
        int yy_lac_status =                                      \
          yy_lac (yyesa, &yyes, &yyes_capacity, yyssp, yytoken); \
        if (yy_lac_status == 2)                                  \
          goto yyexhaustedlab;                                   \
        if (yy_lac_status == 1)                                  \
          goto yyerrlab;                                         \
      }                                                          \
    }                                                            \
} while (0)

/* Discard any previous initial lookahead context because of Event,
   which may be a lookahead change or an invalidation of the currently
   established initial context for the current lookahead.

   The most common example of a lookahead change is a shift.  An example
   of both cases is syntax error recovery.  That is, a syntax error
   occurs when the lookahead is syntactically erroneous for the
   currently established initial context, so error recovery manipulates
   the parser stacks to try to find a new initial context in which the
   current lookahead is syntactically acceptable.  If it fails to find
   such a context, it discards the lookahead.  */
#if YYDEBUG
# define YY_LAC_DISCARD(Event)                                           \
do {                                                                     \
  if (yy_lac_established)                                                \
    {                                                                    \
      if (yydebug)                                                       \
        YYFPRINTF (stderr, "LAC: initial context discarded due to "      \
                   Event "\n");                                          \
      yy_lac_established = 0;                                            \
    }                                                                    \
} while (0)
#else
# define YY_LAC_DISCARD(Event) yy_lac_established = 0
#endif

/* Given the stack whose top is *YYSSP, return 0 iff YYTOKEN can
   eventually (after perhaps some reductions) be shifted, return 1 if
   not, or return 2 if memory is exhausted.  As preconditions and
   postconditions: *YYES_CAPACITY is the allocated size of the array to
   which *YYES points, and either *YYES = YYESA or *YYES points to an
   array allocated with YYSTACK_ALLOC.  yy_lac may overwrite the
   contents of either array, alter *YYES and *YYES_CAPACITY, and free
   any old *YYES other than YYESA.  */
static int
yy_lac (yytype_int16 *yyesa, yytype_int16 **yyes,
        YYSIZE_T *yyes_capacity, yytype_int16 *yyssp, int yytoken)
{
  yytype_int16 *yyes_prev = yyssp;
  yytype_int16 *yyesp = yyes_prev;
  YYDPRINTF ((stderr, "LAC: checking lookahead %s:", yytname[yytoken]));
  if (yytoken == YYUNDEFTOK)
    {
      YYDPRINTF ((stderr, " Always Err\n"));
      return 1;
    }
  while (1)
    {
      int yyrule = yypact[*yyesp];
      if (yypact_value_is_default (yyrule)
          || (yyrule += yytoken) < 0 || YYLAST < yyrule
          || yycheck[yyrule] != yytoken)
        {
          yyrule = yydefact[*yyesp];
          if (yyrule == 0)
            {
              YYDPRINTF ((stderr, " Err\n"));
              return 1;
            }
        }
      else
        {
          yyrule = yytable[yyrule];
          if (yytable_value_is_error (yyrule))
            {
              YYDPRINTF ((stderr, " Err\n"));
              return 1;
            }
          if (0 < yyrule)
            {
              YYDPRINTF ((stderr, " S%d\n", yyrule));
              return 0;
            }
          yyrule = -yyrule;
        }
      {
        YYSIZE_T yylen = yyr2[yyrule];
        YYDPRINTF ((stderr, " R%d", yyrule - 1));
        if (yyesp != yyes_prev)
          {
            YYSIZE_T yysize = (YYSIZE_T) (yyesp - *yyes + 1);
            if (yylen < yysize)
              {
                yyesp -= yylen;
                yylen = 0;
              }
            else
              {
                yylen -= yysize;
                yyesp = yyes_prev;
              }
          }
        if (yylen)
          yyesp = yyes_prev -= yylen;
      }
      {
        yytype_int16 yystate;
        {
          const int yylhs = yyr1[yyrule] - YYNTOKENS;
          const int yyi = yypgoto[yylhs] + *yyesp;
          yystate = ((yytype_int16)
                     (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyesp
                      ? yytable[yyi]
                      : yydefgoto[yylhs]));
        }
        if (yyesp == yyes_prev)
          {
            yyesp = *yyes;
            *yyesp = yystate;
          }
        else
          {
            if (yy_lac_stack_realloc (yyes_capacity, 1,
#if YYDEBUG
                                      " (", ")",
#endif
                                      yyes, yyesa, &yyesp, yyes_prev))
              {
                YYDPRINTF ((stderr, "\n"));
                return 2;
              }
            *++yyesp = yystate;
          }
        YYDPRINTF ((stderr, " G%d", (int) yystate));
      }
    }
}


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

  return (YYSIZE_T) (yystpcpy (yyres, yystr) - yyres);
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.  In order to see if a particular token T is a
   valid looakhead, invoke yy_lac (YYESA, YYES, YYES_CAPACITY, YYSSP, T).

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store or if
   yy_lac returned 2.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyesa, yytype_int16 **yyes,
                YYSIZE_T *yyes_capacity, yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 100 };
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
       In the first two cases, it might appear that the current syntax
       error should have been detected in the previous state when yy_lac
       was invoked.  However, at that time, there might have been a
       different syntax error that discarded a different initial context
       during error recovery, leaving behind the current lookahead.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      YYDPRINTF ((stderr, "Constructing syntax error message\n"));
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          int yyx;

          for (yyx = 0; yyx < YYNTOKENS; ++yyx)
            if (yyx != YYTERROR && yyx != YYUNDEFTOK)
              {
                {
                  int yy_lac_status = yy_lac (yyesa, yyes, yyes_capacity,
                                              yyssp, yyx);
                  if (yy_lac_status == 2)
                    return 2;
                  if (yy_lac_status == 1)
                    continue;
                }
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
# if YYDEBUG
      else if (yydebug)
        YYFPRINTF (stderr, "No expected tokens.\n");
# endif
    }

#define YY_BASE_ERROR_STRING "syntax error, unexpected %s, expecting %s"
#define YY_EXTEND_ERROR_STRING " or %s"
char yyformatbuf[sizeof(YY_BASE_ERROR_STRING) + (sizeof(YY_EXTEND_ERROR_STRING)) * (YYERROR_VERBOSE_ARGS_MAXIMUM - 1)];
strcpy(yyformatbuf, YY_BASE_ERROR_STRING);
yyformat = yyformatbuf + sizeof(YY_BASE_ERROR_STRING) - 1;
{
  int yyi;
  for (yyi = 1; yyi < yycount - 1; yyi++)
  {
    strcpy((char*)yyformat, YY_EXTEND_ERROR_STRING);
    yyformat += sizeof(YY_EXTEND_ERROR_STRING) - 1;
  }
  yyformat = yyformatbuf;
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
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, YYLTYPE *yylocationp, nissa::driver_t *driver)
{
  YYUSE (yyvaluep);
  YYUSE (yylocationp);
  YYUSE (driver);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/*----------.
| yyparse.  |
`----------*/

int
yyparse (nissa::driver_t *driver)
{
/* The lookahead symbol.  */
int yychar;


/* The semantic value of the lookahead symbol.  */
/* Default value used for initialization, for pacifying older GCCs
   or non-GCC compilers.  */
YY_INITIAL_VALUE (static YYSTYPE yyval_default;)
YYSTYPE yylval YY_INITIAL_VALUE (= yyval_default);

/* Location data for the lookahead symbol.  */
static YYLTYPE yyloc_default
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
  = { 1, 1, 1, 1 }
# endif
;
YYLTYPE yylloc = yyloc_default;

    /* Number of syntax errors so far.  */
    int yynerrs;

    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.
       'yyls': related to locations.

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

    /* The location stack.  */
    YYLTYPE yylsa[YYINITDEPTH];
    YYLTYPE *yyls;
    YYLTYPE *yylsp;

    /* The locations where the error started and ended.  */
    YYLTYPE yyerror_range[3];

    YYSIZE_T yystacksize;

    yytype_int16 yyesa[20];
    yytype_int16 *yyes;
    YYSIZE_T yyes_capacity;

  int yy_lac_established = 0;
  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;
  YYLTYPE yyloc;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N), yylsp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yylsp = yyls = yylsa;
  yystacksize = YYINITDEPTH;

  yyes = yyesa;
  yyes_capacity = sizeof yyesa / sizeof *yyes;
  if (YYMAXDEPTH < yyes_capacity)
    yyes_capacity = YYMAXDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  yylsp[0] = yylloc;
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = (yytype_int16) yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = (YYSIZE_T) (yyssp - yyss + 1);

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;
        YYLTYPE *yyls1 = yyls;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yyls1, yysize * sizeof (*yylsp),
                    &yystacksize);
        yyss = yyss1;
        yyvs = yyvs1;
        yyls = yyls1;
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
        YYSTACK_RELOCATE (yyls_alloc, yyls);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;
      yylsp = yyls + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long) yystacksize));

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
      yychar = yylex (&yylval, &yylloc, scanner);
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
    {
      YY_LAC_ESTABLISH;
      goto yydefault;
    }
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      YY_LAC_ESTABLISH;
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
  YY_LAC_DISCARD ("shift");

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END
  *++yylsp = yylloc;
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

  /* Default location. */
  YYLLOC_DEFAULT (yyloc, (yylsp - yylen), yylen);
  yyerror_range[1] = yyloc;
  YY_REDUCE_PRINT (yyn);
  {
    int yychar_backup = yychar;
    switch (yyn)
      {
  case 8:
#line 329 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->tag=(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 2372 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 9:
#line 331 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->theories.push_back(*(yyvsp[0].theory));delete (yyvsp[0].theory);}
#line 2378 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 10:
#line 333 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->add_meson_corr_meas(*(yyvsp[0].meson_corr_meas));delete (yyvsp[0].meson_corr_meas);}
#line 2384 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 11:
#line 334 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->add_nucleon_corr_meas(*(yyvsp[0].nucleon_corr_meas));delete (yyvsp[0].nucleon_corr_meas);}
#line 2390 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 12:
#line 335 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->add_fermionic_putpourri_meas(*(yyvsp[0].fermionic_putpourri_meas));delete (yyvsp[0].fermionic_putpourri_meas);}
#line 2396 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 13:
#line 336 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->add_quark_rendens_meas(*(yyvsp[0].quark_rendens_meas));delete (yyvsp[0].quark_rendens_meas);}
#line 2402 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 14:
#line 337 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->add_chir_zumba_meas(*(yyvsp[0].chir_zumba_meas));delete (yyvsp[0].chir_zumba_meas);}
#line 2408 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 15:
#line 338 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->add_spinpol_meas(*(yyvsp[0].spinpol_meas));delete (yyvsp[0].spinpol_meas);}
#line 2414 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 16:
#line 339 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->add_qed_corr_meas(*(yyvsp[0].qed_corr_meas));delete (yyvsp[0].qed_corr_meas);}
#line 2420 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 17:
#line 340 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->add_magnetization_meas(*(yyvsp[0].magnetization_meas));delete (yyvsp[0].magnetization_meas);}
#line 2426 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 18:
#line 341 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->add_minmax_eigenvalues_meas(*(yyvsp[0].minmax_eigenvalues_meas));delete (yyvsp[0].minmax_eigenvalues_meas);}
#line 2432 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 19:
#line 343 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->add_plaq_pol_meas(*(yyvsp[0].plaq_pol_meas));delete (yyvsp[0].plaq_pol_meas);}
#line 2438 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 20:
#line 344 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->add_top_meas(*(yyvsp[0].top_meas));delete (yyvsp[0].top_meas);}
#line 2444 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 21:
#line 345 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->add_spectr_proj_meas(*(yyvsp[0].spectr_proj_meas));delete (yyvsp[0].spectr_proj_meas);}
#line 2450 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 22:
#line 346 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->add_luppoli_meas(*(yyvsp[0].luppoli_meas));delete (yyvsp[0].luppoli_meas);}
#line 2456 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 23:
#line 347 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->add_watusso_meas(*(yyvsp[0].watusso_meas));delete (yyvsp[0].watusso_meas);}
#line 2462 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 24:
#line 348 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->add_all_rects_meas(*(yyvsp[0].all_rects_meas));delete (yyvsp[0].all_rects_meas);}
#line 2468 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 29:
#line 358 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {FILE *fout=open_file((yyvsp[-1].text)->c_str(),"w");driver->master_fprintf(fout);close_file(fout);delete (yyvsp[-1].text);}
#line 2474 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 30:
#line 359 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {FILE *fout=open_file((yyvsp[-1].text)->c_str(),"w");driver->master_fprintf(fout,true);close_file(fout);delete (yyvsp[-1].text);}
#line 2480 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 31:
#line 360 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {master_printf((yyvsp[-1].text)->c_str());delete (yyvsp[-1].text);}
#line 2486 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 32:
#line 365 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.gauge_action_name)=WILSON_GAUGE_ACTION;}
#line 2492 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 33:
#line 366 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.gauge_action_name)=TLSYM_GAUGE_ACTION;}
#line 2498 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 34:
#line 367 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.gauge_action_name)=IWASAKI_GAUGE_ACTION;}
#line 2504 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 35:
#line 369 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_numb)=(yyvsp[0].int_numb);}
#line 2510 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 36:
#line 371 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_numb)=(yyvsp[0].double_numb);}
#line 2516 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 37:
#line 373 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[0].int_numb);}
#line 2522 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 38:
#line 375 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_numb)=(yyvsp[0].double_numb);}
#line 2528 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 39:
#line 377 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[0].int_numb);}
#line 2534 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 40:
#line 379 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[0].int_numb);}
#line 2540 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 41:
#line 381 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[0].int_numb);}
#line 2546 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 42:
#line 383 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.rnd_type)=(yyvsp[0].rnd_type);}
#line 2552 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 43:
#line 385 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[0].int_numb);}
#line 2558 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 44:
#line 387 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[0].int_numb);}
#line 2564 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 45:
#line 389 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_numb)=(yyvsp[0].double_numb);}
#line 2570 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 46:
#line 391 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[0].int_numb);}
#line 2576 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 47:
#line 393 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[0].int_numb);}
#line 2582 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 48:
#line 395 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[0].int_numb);}
#line 2588 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 49:
#line 397 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.text)=(yyvsp[0].text);}
#line 2594 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 51:
#line 400 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->T=(yyvsp[0].int_numb);}
#line 2600 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 52:
#line 401 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->LX=driver->LY=driver->LZ=(yyvsp[0].int_numb);}
#line 2606 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 53:
#line 402 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->LX=(yyvsp[0].int_numb);}
#line 2612 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 54:
#line 403 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->LY=(yyvsp[0].int_numb);}
#line 2618 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 55:
#line 404 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->LZ=(yyvsp[0].int_numb);}
#line 2624 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 56:
#line 409 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.theory)=new theory_pars_t;}
#line 2630 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 57:
#line 410 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.theory)->stout_pars=(*(yyvsp[0].stout_pars));delete (yyvsp[0].stout_pars);}
#line 2636 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 58:
#line 411 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyvsp[0].topotential_pars)->init();(yyval.theory)->topotential_pars=(*(yyvsp[0].topotential_pars));delete (yyvsp[0].topotential_pars);}
#line 2642 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 59:
#line 412 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.theory)->quarks.push_back(*(yyvsp[0].quark));delete (yyvsp[0].quark);}
#line 2648 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 60:
#line 413 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.theory)->em_field_pars=(*(yyvsp[0].em_field_pars));delete (yyvsp[0].em_field_pars);}
#line 2654 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 61:
#line 414 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.theory)->beta=(yyvsp[0].double_numb);}
#line 2660 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 62:
#line 415 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.theory)->gauge_action_name=(yyvsp[0].gauge_action_name);}
#line 2666 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 63:
#line 419 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.em_field_pars)=new em_field_pars_t;(yyval.em_field_pars)->flag=1;}
#line 2672 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 64:
#line 420 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.em_field_pars)->E[(yyvsp[-2].int_numb)]=(yyvsp[0].double_numb);}
#line 2678 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 65:
#line 421 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.em_field_pars)->B[(yyvsp[-2].int_numb)]=(yyvsp[0].double_numb);}
#line 2684 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 66:
#line 426 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {}
#line 2690 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 67:
#line 427 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->walltime=(yyvsp[0].int_numb);}
#line 2696 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 68:
#line 428 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->seed=(yyvsp[0].int_numb);}
#line 2702 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 69:
#line 433 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.topotential_pars)=new topotential_pars_t;(yyval.topotential_pars)->flag=2;}
#line 2708 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 70:
#line 434 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.topotential_pars)=new topotential_pars_t;(yyval.topotential_pars)->flag=1;}
#line 2714 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 71:
#line 435 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.topotential_pars)=new topotential_pars_t;(yyval.topotential_pars)->flag=0;}
#line 2720 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 72:
#line 436 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.topotential_pars)->stout_pars=(*(yyvsp[0].stout_pars));delete (yyvsp[0].stout_pars);}
#line 2726 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 73:
#line 437 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.topotential_pars)->after=(yyvsp[0].int_numb);}
#line 2732 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 74:
#line 438 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.topotential_pars)->each=(yyvsp[0].double_numb);}
#line 2738 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 75:
#line 439 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.topotential_pars)->theta=(yyvsp[0].double_numb);}
#line 2744 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 76:
#line 440 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.topotential_pars)->coeff=(yyvsp[0].double_numb);}
#line 2750 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 77:
#line 441 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.topotential_pars)->width=(yyvsp[0].double_numb);}
#line 2756 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 78:
#line 442 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.topotential_pars)->barr=(yyvsp[0].double_numb);}
#line 2762 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 79:
#line 443 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.topotential_pars)->force_out=(yyvsp[0].double_numb);}
#line 2768 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 80:
#line 444 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.topotential_pars)->well_tempering=(yyvsp[0].double_numb);}
#line 2774 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 81:
#line 445 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.topotential_pars)->bend=(yyvsp[0].double_numb);}
#line 2780 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 82:
#line 450 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.smooth_pars)=new smooth_pars_t;(yyval.smooth_pars)->method=smooth_pars_t::STOUT;(yyval.smooth_pars)->stout=(*(yyvsp[0].stout_pars));delete (yyvsp[0].stout_pars);}
#line 2786 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 83:
#line 451 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.smooth_pars)=new smooth_pars_t;(yyval.smooth_pars)->method=smooth_pars_t::COOLING;(yyval.smooth_pars)->cool=(*(yyvsp[0].cool_pars));delete (yyvsp[0].cool_pars);}
#line 2792 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 84:
#line 452 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.smooth_pars)=new smooth_pars_t;(yyval.smooth_pars)->method=smooth_pars_t::WFLOW;(yyval.smooth_pars)->Wflow=(*(yyvsp[0].Wflow_pars));delete (yyvsp[0].Wflow_pars);}
#line 2798 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 85:
#line 453 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.smooth_pars)=new smooth_pars_t;(yyval.smooth_pars)->method=smooth_pars_t::APE;(yyval.smooth_pars)->ape=(*(yyvsp[0].ape_pars));delete (yyvsp[0].ape_pars);}
#line 2804 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 86:
#line 454 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.smooth_pars)=new smooth_pars_t;(yyval.smooth_pars)->method=smooth_pars_t::HYP;(yyval.smooth_pars)->hyp=(*(yyvsp[0].hyp_pars));delete (yyvsp[0].hyp_pars);}
#line 2810 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 87:
#line 455 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyvsp[-1].smooth_pars)->meas_each_nsmooth=(yyvsp[0].double_numb);}
#line 2816 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 88:
#line 456 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.smooth_pars)->space_or_time=smooth_pars_t::SPACE;}
#line 2822 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 89:
#line 457 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.smooth_pars)->space_or_time=smooth_pars_t::TIME;}
#line 2828 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 90:
#line 458 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.smooth_pars)->space_or_time=smooth_pars_t::SPACETIME;}
#line 2834 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 91:
#line 463 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.stout_pars)=new stout_pars_t;}
#line 2840 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 92:
#line 464 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.stout_pars)=(yyvsp[-1].stout_pars);(yyval.stout_pars)->nlevels=(yyvsp[0].int_numb);}
#line 2846 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 93:
#line 465 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.stout_pars)=(yyvsp[-1].stout_pars);(yyval.stout_pars)->rho=(yyvsp[0].double_numb);}
#line 2852 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 94:
#line 468 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[0].int_numb);}
#line 2858 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 95:
#line 469 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_numb)=(yyvsp[0].double_numb);}
#line 2864 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 96:
#line 473 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.cool_pars)=new cool_pars_t;}
#line 2870 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 97:
#line 474 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.cool_pars)=(yyvsp[-1].cool_pars);(yyval.cool_pars)->nsteps=(yyvsp[0].int_numb);}
#line 2876 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 98:
#line 475 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.cool_pars)=(yyvsp[-1].cool_pars);(yyval.cool_pars)->gauge_action=(yyvsp[0].gauge_action_name);}
#line 2882 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 99:
#line 478 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[0].int_numb);}
#line 2888 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 100:
#line 482 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.ape_pars)=new ape_pars_t;}
#line 2894 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 101:
#line 483 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.ape_pars)=(yyvsp[-1].ape_pars);(yyval.ape_pars)->nlevels=(yyvsp[0].int_numb);}
#line 2900 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 102:
#line 484 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.ape_pars)=(yyvsp[-1].ape_pars);(yyval.ape_pars)->alpha=(yyvsp[0].double_numb);}
#line 2906 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 103:
#line 487 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_numb)=(yyvsp[0].double_numb);}
#line 2912 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 104:
#line 491 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.hyp_pars)=new hyp_pars_t;}
#line 2918 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 105:
#line 492 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.hyp_pars)=(yyvsp[-1].hyp_pars);(yyval.hyp_pars)->nlevels=(yyvsp[0].int_numb);}
#line 2924 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 106:
#line 493 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {
	    (yyval.hyp_pars)=(yyvsp[-1].hyp_pars);
	    if((yyvsp[0].double_list)->size()!=3) crash("hyp alpha needs an exactly 3 long list of int");
	    (yyval.hyp_pars)->alpha0=(*(yyvsp[0].double_list))[0];
	    (yyval.hyp_pars)->alpha1=(*(yyvsp[0].double_list))[1];
	    (yyval.hyp_pars)->alpha2=(*(yyvsp[0].double_list))[2];
	    delete (yyvsp[0].double_list);
	    }
#line 2937 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 107:
#line 503 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_list)=(yyvsp[0].double_list);}
#line 2943 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 108:
#line 507 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.Wflow_pars)=new Wflow_pars_t;}
#line 2949 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 109:
#line 508 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.Wflow_pars)=(yyvsp[-1].Wflow_pars);(yyval.Wflow_pars)->nflows=(yyvsp[0].int_numb);}
#line 2955 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 110:
#line 509 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.Wflow_pars)=(yyvsp[-1].Wflow_pars);(yyval.Wflow_pars)->nrecu=(yyvsp[0].int_numb);}
#line 2961 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 111:
#line 510 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.Wflow_pars)=(yyvsp[-1].Wflow_pars);(yyval.Wflow_pars)->dt=(yyvsp[0].double_numb);}
#line 2967 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 112:
#line 513 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[0].int_numb);}
#line 2973 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 113:
#line 514 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[0].int_numb);}
#line 2979 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 114:
#line 515 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_numb)=(yyvsp[0].double_numb);}
#line 2985 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 115:
#line 520 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {
	 (yyval.quark)=new quark_content_t();
	 (yyval.quark)->name=(*(yyvsp[0].text));
	 delete (yyvsp[0].text);
     }
#line 2995 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 116:
#line 525 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark)->deg=(yyvsp[0].int_numb);}
#line 3001 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 117:
#line 526 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark)->discretiz=nissa::ferm_discretiz::ROOT_STAG;}
#line 3007 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 118:
#line 527 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark)->discretiz=nissa::ferm_discretiz::ROOT_TM_CLOV;}
#line 3013 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 119:
#line 528 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark)->discretiz=nissa::ferm_discretiz::OVERLAP;}
#line 3019 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 120:
#line 529 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark)->mass=(yyvsp[0].double_numb);}
#line 3025 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 121:
#line 530 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark)->mass_overlap=(yyvsp[0].double_numb);}
#line 3031 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 122:
#line 531 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark)->kappa=(yyvsp[0].double_numb);}
#line 3037 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 123:
#line 532 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark)->cSW=(yyvsp[0].double_numb);}
#line 3043 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 124:
#line 533 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark)->re_pot=(yyvsp[0].double_numb);}
#line 3049 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 125:
#line 534 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark)->im_pot=(yyvsp[0].double_numb);}
#line 3055 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 126:
#line 535 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark)->charge=(yyvsp[0].double_numb);}
#line 3061 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 127:
#line 540 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.nucleon_corr_meas)=new nucleon_corr_meas_pars_t();}
#line 3067 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 128:
#line 541 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.nucleon_corr_meas)->each=(yyvsp[0].double_numb);}
#line 3073 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 129:
#line 542 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.nucleon_corr_meas)->after=(yyvsp[0].int_numb);}
#line 3079 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 130:
#line 543 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.nucleon_corr_meas)->path=(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 3085 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 131:
#line 544 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.nucleon_corr_meas)->residue=(yyvsp[0].double_numb);}
#line 3091 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 132:
#line 545 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.nucleon_corr_meas)->itheory=(yyvsp[0].int_numb);}
#line 3097 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 133:
#line 546 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.nucleon_corr_meas)->ncopies=(yyvsp[0].int_numb);}
#line 3103 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 134:
#line 547 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.nucleon_corr_meas)->rnd_type=(yyvsp[0].rnd_type);}
#line 3109 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 135:
#line 548 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.nucleon_corr_meas)->nhits=(yyvsp[0].int_numb);}
#line 3115 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 136:
#line 553 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.meson_corr_meas)=new meson_corr_meas_pars_t();}
#line 3121 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 137:
#line 554 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.meson_corr_meas)->each=(yyvsp[0].double_numb);}
#line 3127 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 138:
#line 555 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.meson_corr_meas)->after=(yyvsp[0].int_numb);}
#line 3133 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 139:
#line 556 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.meson_corr_meas)->path=(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 3139 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 140:
#line 557 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.meson_corr_meas)->residue=(yyvsp[0].double_numb);}
#line 3145 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 141:
#line 558 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.meson_corr_meas)->itheory=(yyvsp[0].int_numb);}
#line 3151 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 142:
#line 559 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.meson_corr_meas)->ncopies=(yyvsp[0].int_numb);}
#line 3157 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 143:
#line 560 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.meson_corr_meas)->rnd_type=(yyvsp[0].rnd_type);}
#line 3163 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 144:
#line 561 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.meson_corr_meas)->nhits=(yyvsp[0].int_numb);}
#line 3169 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 145:
#line 562 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {
		   if((*(yyvsp[0].int_pair_list))[0]!=std::make_pair(15,15)) crash("first entry of the meson list must be 15,15, instead it is: %d %d",(*(yyvsp[0].int_pair_list))[0].first,(*(yyvsp[0].int_pair_list))[0].second);
		   (yyval.meson_corr_meas)->mesons=(*(yyvsp[0].int_pair_list));delete (yyvsp[0].int_pair_list);}
#line 3177 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 146:
#line 567 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {
    for(size_t i=0;i<(yyvsp[0].int_pair_list)->size();i++)
	{
	    int f=(*(yyvsp[0].int_pair_list))[i].first;
	    int s=(*(yyvsp[0].int_pair_list))[i].second;
	    if((f<0)||(f>15)) crash("first part of entry %d is %d, should be in the range [0,15]",f);
	    if((s<0)||(s>15)) crash("second part of entry %d is %d, should be in the range [0,15]",s);
	}
	(yyval.int_pair_list)=(yyvsp[0].int_pair_list);}
#line 3191 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 147:
#line 581 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.fermionic_putpourri_meas)=new fermionic_putpourri_meas_pars_t();}
#line 3197 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 148:
#line 582 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.fermionic_putpourri_meas)->each=(yyvsp[0].double_numb);}
#line 3203 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 149:
#line 583 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.fermionic_putpourri_meas)->after=(yyvsp[0].int_numb);}
#line 3209 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 150:
#line 584 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.fermionic_putpourri_meas)->path=(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 3215 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 151:
#line 585 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.fermionic_putpourri_meas)->residue=(yyvsp[0].double_numb);}
#line 3221 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 152:
#line 586 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.fermionic_putpourri_meas)->compute_susc=(yyvsp[0].int_numb);}
#line 3227 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 153:
#line 587 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.fermionic_putpourri_meas)->itheory=(yyvsp[0].int_numb);}
#line 3233 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 154:
#line 588 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.fermionic_putpourri_meas)->ncopies=(yyvsp[0].int_numb);}
#line 3239 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 155:
#line 589 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.fermionic_putpourri_meas)->rnd_type=(yyvsp[0].rnd_type);}
#line 3245 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 156:
#line 590 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.fermionic_putpourri_meas)->nhits=(yyvsp[0].int_numb);}
#line 3251 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 157:
#line 595 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark_rendens_meas)=new quark_rendens_meas_pars_t();}
#line 3257 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 158:
#line 596 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark_rendens_meas)->each=(yyvsp[0].double_numb);}
#line 3263 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 159:
#line 597 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark_rendens_meas)->after=(yyvsp[0].int_numb);}
#line 3269 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 160:
#line 598 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark_rendens_meas)->path=(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 3275 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 161:
#line 599 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark_rendens_meas)->residue=(yyvsp[0].double_numb);}
#line 3281 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 162:
#line 600 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark_rendens_meas)->itheory=(yyvsp[0].int_numb);}
#line 3287 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 163:
#line 601 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark_rendens_meas)->ncopies=(yyvsp[0].int_numb);}
#line 3293 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 164:
#line 602 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark_rendens_meas)->rnd_type=(yyvsp[0].rnd_type);}
#line 3299 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 165:
#line 603 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark_rendens_meas)->nhits=(yyvsp[0].int_numb);}
#line 3305 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 166:
#line 604 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.quark_rendens_meas)->max_order=(yyvsp[0].int_numb);}
#line 3311 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 167:
#line 609 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.chir_zumba_meas)=new chir_zumba_meas_pars_t();}
#line 3317 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 168:
#line 610 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.chir_zumba_meas)->each=(yyvsp[0].double_numb);}
#line 3323 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 169:
#line 611 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.chir_zumba_meas)->after=(yyvsp[0].int_numb);}
#line 3329 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 170:
#line 612 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.chir_zumba_meas)->path=(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 3335 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 171:
#line 613 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.chir_zumba_meas)->residue=(yyvsp[0].double_numb);}
#line 3341 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 172:
#line 614 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.chir_zumba_meas)->itheory=(yyvsp[0].int_numb);}
#line 3347 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 173:
#line 615 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.chir_zumba_meas)->ncopies=(yyvsp[0].int_numb);}
#line 3353 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 174:
#line 616 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.chir_zumba_meas)->rnd_type=(yyvsp[0].rnd_type);}
#line 3359 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 175:
#line 617 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.chir_zumba_meas)->nhits=(yyvsp[0].int_numb);}
#line 3365 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 176:
#line 618 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.chir_zumba_meas)->max_order=(yyvsp[0].int_numb);}
#line 3371 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 177:
#line 623 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spinpol_meas)=new spinpol_meas_pars_t();}
#line 3377 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 178:
#line 624 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spinpol_meas)->each=(yyvsp[0].double_numb);}
#line 3383 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 179:
#line 625 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spinpol_meas)->after=(yyvsp[0].int_numb);}
#line 3389 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 180:
#line 626 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spinpol_meas)->path=(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 3395 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 181:
#line 627 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spinpol_meas)->residue=(yyvsp[0].double_numb);}
#line 3401 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 182:
#line 628 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spinpol_meas)->itheory=(yyvsp[0].int_numb);}
#line 3407 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 183:
#line 629 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spinpol_meas)->ncopies=(yyvsp[0].int_numb);}
#line 3413 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 184:
#line 630 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spinpol_meas)->rnd_type=(yyvsp[0].rnd_type);}
#line 3419 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 185:
#line 631 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spinpol_meas)->nhits=(yyvsp[0].int_numb);}
#line 3425 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 186:
#line 632 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spinpol_meas)->operators=*(yyvsp[0].int_pair_list);delete (yyvsp[0].int_pair_list);}
#line 3431 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 187:
#line 633 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spinpol_meas)->use_ferm_conf_for_gluons=(yyvsp[0].int_numb);}
#line 3437 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 188:
#line 634 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spinpol_meas)->use_adjoint_flow=(yyvsp[0].int_numb);}
#line 3443 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 189:
#line 635 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spinpol_meas)->smooth_pars=(*(yyvsp[0].smooth_pars));delete (yyvsp[0].smooth_pars);}
#line 3449 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 190:
#line 638 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[0].int_numb);}
#line 3455 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 191:
#line 639 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[0].int_numb);}
#line 3461 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 192:
#line 643 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.qed_corr_meas)=new qed_corr_meas_pars_t();}
#line 3467 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 193:
#line 644 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.qed_corr_meas)->each=(yyvsp[0].double_numb);}
#line 3473 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 194:
#line 645 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.qed_corr_meas)->after=(yyvsp[0].int_numb);}
#line 3479 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 195:
#line 646 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.qed_corr_meas)->path=(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 3485 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 196:
#line 647 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.qed_corr_meas)->residue=(yyvsp[0].double_numb);}
#line 3491 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 197:
#line 648 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.qed_corr_meas)->itheory=(yyvsp[0].int_numb);}
#line 3497 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 198:
#line 649 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.qed_corr_meas)->ncopies=(yyvsp[0].int_numb);}
#line 3503 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 199:
#line 650 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.qed_corr_meas)->rnd_type=(yyvsp[0].rnd_type);}
#line 3509 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 200:
#line 651 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.qed_corr_meas)->nhits=(yyvsp[0].int_numb);}
#line 3515 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 201:
#line 656 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.magnetization_meas)=new magnetization_meas_pars_t();}
#line 3521 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 202:
#line 657 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.magnetization_meas)->each=(yyvsp[0].double_numb);}
#line 3527 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 203:
#line 658 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.magnetization_meas)->after=(yyvsp[0].int_numb);}
#line 3533 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 204:
#line 659 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.magnetization_meas)->path=(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 3539 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 205:
#line 660 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.magnetization_meas)->residue=(yyvsp[0].double_numb);}
#line 3545 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 206:
#line 661 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.magnetization_meas)->itheory=(yyvsp[0].int_numb);}
#line 3551 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 207:
#line 662 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.magnetization_meas)->ncopies=(yyvsp[0].int_numb);}
#line 3557 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 208:
#line 663 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.magnetization_meas)->rnd_type=(yyvsp[0].rnd_type);}
#line 3563 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 209:
#line 664 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.magnetization_meas)->nhits=(yyvsp[0].int_numb);}
#line 3569 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 210:
#line 669 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.minmax_eigenvalues_meas)=new minmax_eigenvalues_meas_pars_t();}
#line 3575 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 211:
#line 670 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.minmax_eigenvalues_meas)->each=(yyvsp[0].double_numb);}
#line 3581 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 212:
#line 671 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.minmax_eigenvalues_meas)->after=(yyvsp[0].int_numb);}
#line 3587 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 213:
#line 672 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.minmax_eigenvalues_meas)->path=(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 3593 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 214:
#line 673 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.minmax_eigenvalues_meas)->residue=(yyvsp[0].double_numb);}
#line 3599 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 215:
#line 674 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.minmax_eigenvalues_meas)->itheory=(yyvsp[0].int_numb);}
#line 3605 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 216:
#line 675 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.minmax_eigenvalues_meas)->ncopies=(yyvsp[0].int_numb);}
#line 3611 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 217:
#line 676 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.minmax_eigenvalues_meas)->rnd_type=(yyvsp[0].rnd_type);}
#line 3617 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 218:
#line 677 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.minmax_eigenvalues_meas)->nhits=(yyvsp[0].int_numb);}
#line 3623 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 219:
#line 678 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.minmax_eigenvalues_meas)->neigs=(yyvsp[0].int_numb);}
#line 3629 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 220:
#line 679 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.minmax_eigenvalues_meas)->wspace_size=(yyvsp[0].int_numb);}
#line 3635 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 221:
#line 680 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.minmax_eigenvalues_meas)->min_max=(yyvsp[0].int_numb);}
#line 3641 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 222:
#line 685 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->run_mode=driver_t::EVOLUTION_MODE;}
#line 3647 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 223:
#line 686 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->evol_pars.ntraj_tot=(yyvsp[0].int_numb);}
#line 3653 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 224:
#line 687 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->evol_pars.id_sea_theory=(yyvsp[0].int_numb);}
#line 3659 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 225:
#line 688 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->evol_pars.skip_mtest_ntraj=(yyvsp[0].int_numb);}
#line 3665 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 226:
#line 689 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->evol_pars.traj_length=(yyvsp[0].double_numb);}
#line 3671 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 227:
#line 690 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->evol_pars.pf_action_residue=(yyvsp[0].double_numb);}
#line 3677 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 228:
#line 691 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->evol_pars.md_residue=(yyvsp[0].double_numb);}
#line 3683 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 229:
#line 692 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->evol_pars.nmd_steps=(yyvsp[0].int_numb);}
#line 3689 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 230:
#line 693 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->evol_pars.ngauge_substeps=(yyvsp[0].int_numb);}
#line 3695 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 231:
#line 694 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->evol_pars.npseudo_fs=(*(yyvsp[0].int_list));delete (yyvsp[0].int_list);}
#line 3701 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 232:
#line 699 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->run_mode=driver_t::EVOLUTION_MODE;}
#line 3707 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 233:
#line 700 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->conf_pars.path=(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 3713 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 234:
#line 702 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {
                        driver->conf_pars.store_path=(*(yyvsp[0].text));
                        char test1[128];
                        char test2[128];
                        snprintf(test1,128,(yyvsp[0].text)->c_str(),100);
                        snprintf(test2,128,(yyvsp[0].text)->c_str(),101);
                        if(!strcmp(test1,test2)) crash("bad template \"%s\" for store_path",(yyvsp[0].text)->c_str());
			delete (yyvsp[0].text);
		    }
#line 3727 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 235:
#line 711 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->conf_pars.store_each=(yyvsp[0].int_numb);}
#line 3733 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 236:
#line 712 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->conf_pars.store_running=(yyvsp[0].int_numb);}
#line 3739 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 237:
#line 713 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->conf_pars.start_cond=HOT_START_COND;}
#line 3745 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 238:
#line 714 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {driver->conf_pars.start_cond=COLD_START_COND;}
#line 3751 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 239:
#line 719 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spectr_proj_meas)=new spectr_proj_meas_pars_t();}
#line 3757 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 240:
#line 720 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spectr_proj_meas)->each=(yyvsp[0].double_numb);}
#line 3763 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 241:
#line 721 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spectr_proj_meas)->after=(yyvsp[0].int_numb);}
#line 3769 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 242:
#line 722 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spectr_proj_meas)->path=(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 3775 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 243:
#line 723 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spectr_proj_meas)->residue=(yyvsp[0].double_numb);}
#line 3781 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 244:
#line 724 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spectr_proj_meas)->itheory=(yyvsp[0].int_numb);}
#line 3787 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 245:
#line 725 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spectr_proj_meas)->ncopies=(yyvsp[0].int_numb);}
#line 3793 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 246:
#line 726 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spectr_proj_meas)->rnd_type=(yyvsp[0].rnd_type);}
#line 3799 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 247:
#line 727 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spectr_proj_meas)->nhits=(yyvsp[0].int_numb);}
#line 3805 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 248:
#line 728 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spectr_proj_meas)->neigs=(yyvsp[0].int_numb);}
#line 3811 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 249:
#line 729 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spectr_proj_meas)->eig_precision=(yyvsp[0].double_numb);}
#line 3817 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 250:
#line 730 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.spectr_proj_meas)->wspace_size=(yyvsp[0].int_numb);}
#line 3823 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 251:
#line 736 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {
    driver->run_mode=driver_t::ANALYSIS_MODE;
    for(auto& it : *(yyvsp[0].text_list))
	{
	    glob_t globbuf;
	    if(glob(it.c_str(),0,NULL,&globbuf))
		crash("Unable to find patterm %s for conf",it.c_str());
	    else
		for(int j=0;j<(int)globbuf.gl_pathc;j++)
		    driver->an_conf_list.push_back(globbuf.gl_pathv[j]);
	    globfree(&globbuf);
	 }
     delete (yyvsp[0].text_list);
 }
#line 3842 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 252:
#line 754 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.plaq_pol_meas)=new gauge_obs_meas_pars_t();}
#line 3848 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 253:
#line 755 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.plaq_pol_meas)->each=(yyvsp[0].double_numb);}
#line 3854 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 254:
#line 756 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.plaq_pol_meas)->after=(yyvsp[0].int_numb);}
#line 3860 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 255:
#line 757 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.plaq_pol_meas)->path=(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 3866 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 256:
#line 758 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.plaq_pol_meas)->use_smooth=(yyvsp[0].int_numb);}
#line 3872 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 257:
#line 759 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.plaq_pol_meas)->meas_plaq=(yyvsp[0].int_numb);}
#line 3878 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 258:
#line 760 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.plaq_pol_meas)->meas_energy=(yyvsp[0].int_numb);}
#line 3884 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 259:
#line 761 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.plaq_pol_meas)->meas_poly=(yyvsp[0].int_numb);}
#line 3890 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 260:
#line 762 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.plaq_pol_meas)->smooth_pars=(*(yyvsp[0].smooth_pars));delete (yyvsp[0].smooth_pars);}
#line 3896 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 261:
#line 767 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.top_meas)=new top_meas_pars_t();}
#line 3902 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 262:
#line 768 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.top_meas)->each=(yyvsp[0].double_numb);}
#line 3908 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 263:
#line 769 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.top_meas)->after=(yyvsp[0].int_numb);}
#line 3914 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 264:
#line 770 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.top_meas)->path=(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 3920 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 265:
#line 771 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.top_meas)->smooth_pars=(*(yyvsp[0].smooth_pars));delete (yyvsp[0].smooth_pars);}
#line 3926 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 266:
#line 772 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.top_meas)->meas_corr=(yyvsp[0].int_numb);}
#line 3932 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 267:
#line 773 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.top_meas)->corr_path=(*(yyvsp[0].text));}
#line 3938 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 268:
#line 778 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.luppoli_meas)=new poly_corr_meas_pars_t();}
#line 3944 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 269:
#line 779 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.luppoli_meas)->each=(yyvsp[0].double_numb);}
#line 3950 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 270:
#line 780 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.luppoli_meas)->after=(yyvsp[0].int_numb);}
#line 3956 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 271:
#line 781 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.luppoli_meas)->path=(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 3962 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 272:
#line 786 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.watusso_meas)=new watusso_meas_pars_t();}
#line 3968 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 273:
#line 787 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.watusso_meas)->each=(yyvsp[0].double_numb);}
#line 3974 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 274:
#line 788 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.watusso_meas)->after=(yyvsp[0].int_numb);}
#line 3980 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 275:
#line 789 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.watusso_meas)->spat_smear_pars=(*(yyvsp[0].smooth_pars));delete (yyvsp[0].smooth_pars);}
#line 3986 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 276:
#line 790 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.watusso_meas)->temp_smear_pars=(*(yyvsp[0].smooth_pars));delete (yyvsp[0].smooth_pars);}
#line 3992 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 277:
#line 791 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.watusso_meas)->path=(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 3998 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 278:
#line 796 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.all_rects_meas)=new all_rects_meas_pars_t();}
#line 4004 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 279:
#line 797 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.all_rects_meas)->each=(yyvsp[0].double_numb);}
#line 4010 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 280:
#line 798 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.all_rects_meas)->after=(yyvsp[0].int_numb);}
#line 4016 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 281:
#line 799 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.all_rects_meas)->Dmin=(yyvsp[0].int_numb);}
#line 4022 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 282:
#line 800 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.all_rects_meas)->Dmax=(yyvsp[0].int_numb);}
#line 4028 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 283:
#line 801 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.all_rects_meas)->Tmin=(yyvsp[0].int_numb);}
#line 4034 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 284:
#line 802 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.all_rects_meas)->Tmax=(yyvsp[0].int_numb);}
#line 4040 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 285:
#line 803 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.all_rects_meas)->spat_smear_pars=(*(yyvsp[0].smooth_pars));delete (yyvsp[0].smooth_pars);}
#line 4046 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 286:
#line 804 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.all_rects_meas)->temp_smear_pars=(*(yyvsp[0].smooth_pars));delete (yyvsp[0].smooth_pars);}
#line 4052 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 287:
#line 805 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.all_rects_meas)->path=(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 4058 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 288:
#line 810 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_list)=(yyvsp[-1].int_list);}
#line 4064 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 289:
#line 813 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_list)=new std::vector<int>;(yyval.int_list)->push_back((yyvsp[0].int_numb));}
#line 4070 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 290:
#line 814 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_list)->push_back((yyvsp[0].int_numb));}
#line 4076 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 291:
#line 819 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_list)=(yyvsp[-1].double_list);}
#line 4082 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 292:
#line 822 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_list)=new std::vector<double>;(yyval.double_list)->push_back((yyvsp[0].double_numb));}
#line 4088 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 293:
#line 823 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_list)->push_back((yyvsp[0].double_numb));}
#line 4094 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 294:
#line 828 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.text_list)=(yyvsp[-1].text_list);}
#line 4100 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 295:
#line 831 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.text_list)=new std::vector<std::string>;(yyval.text_list)->push_back(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 4106 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 296:
#line 832 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.text_list)->push_back(*(yyvsp[0].text));delete (yyvsp[0].text);}
#line 4112 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 297:
#line 837 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_pair_list)=(yyvsp[-1].int_pair_list);}
#line 4118 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 298:
#line 840 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_pair_list)=new std::vector<std::pair<int,int> >;(yyval.int_pair_list)->push_back(*(yyvsp[0].int_pair));delete (yyvsp[0].int_pair);}
#line 4124 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 299:
#line 841 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_pair_list)->push_back(*(yyvsp[0].int_pair));delete (yyvsp[0].int_pair);}
#line 4130 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 300:
#line 846 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_pair)=new std::pair<int,int>; (*(yyval.int_pair))=std::make_pair((yyvsp[-3].int_numb),(yyvsp[-1].int_numb));}
#line 4136 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 301:
#line 852 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.text)=(yyvsp[0].text);}
#line 4142 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 302:
#line 853 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.text)=new std::string(*(yyvsp[-2].text));(*(yyval.text))+=(*(yyvsp[0].text));delete (yyvsp[-2].text);delete (yyvsp[0].text);}
#line 4148 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 303:
#line 857 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_numb)=(yyvsp[0].double_numb);}
#line 4154 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 304:
#line 858 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_numb)=(yyvsp[0].int_numb);}
#line 4160 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 305:
#line 859 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_numb)=(yyvsp[0].double_numb);}
#line 4166 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 306:
#line 860 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_numb)=-(yyvsp[0].double_numb);}
#line 4172 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 307:
#line 861 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_numb)=(yyvsp[-2].double_numb)+(yyvsp[0].double_numb);}
#line 4178 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 308:
#line 862 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_numb)=(yyvsp[-2].double_numb)-(yyvsp[0].double_numb);}
#line 4184 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 309:
#line 863 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_numb)=(yyvsp[-2].double_numb)*(yyvsp[0].double_numb);}
#line 4190 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 310:
#line 864 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_numb)=(yyvsp[-2].double_numb)/(yyvsp[0].double_numb);}
#line 4196 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 311:
#line 865 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_numb)=pow((yyvsp[-2].double_numb),(yyvsp[0].double_numb));}
#line 4202 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 312:
#line 866 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.double_numb)=(yyvsp[-1].double_numb);}
#line 4208 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 313:
#line 870 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[0].int_numb);}
#line 4214 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 314:
#line 871 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[0].int_numb);}
#line 4220 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 315:
#line 872 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=-(yyvsp[0].int_numb);}
#line 4226 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 316:
#line 873 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[-2].int_numb)+(yyvsp[0].int_numb);}
#line 4232 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 317:
#line 874 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[-2].int_numb)-(yyvsp[0].int_numb);}
#line 4238 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 318:
#line 875 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[-1].int_numb)*(yyvsp[0].int_numb);}
#line 4244 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 319:
#line 876 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[-2].int_numb)*(yyvsp[0].int_numb);}
#line 4250 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 320:
#line 877 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[-2].int_numb)/(yyvsp[0].int_numb);}
#line 4256 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 321:
#line 878 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(int)pow((yyvsp[-2].int_numb),(yyvsp[0].int_numb));}
#line 4262 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;

  case 322:
#line 879 "../../projects/generate_confs/parser.ypp" /* yacc.c:1645  */
    {(yyval.int_numb)=(yyvsp[-1].int_numb);}
#line 4268 "generate_confs/parser.cpp" /* yacc.c:1645  */
    break;


#line 4272 "generate_confs/parser.cpp" /* yacc.c:1645  */
        default: break;
      }
    if (yychar_backup != yychar)
      YY_LAC_DISCARD ("yychar change");
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
  *++yylsp = yyloc;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */
  {
    const int yylhs = yyr1[yyn] - YYNTOKENS;
    const int yyi = yypgoto[yylhs] + *yyssp;
    yystate = (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyssp
               ? yytable[yyi]
               : yydefgoto[yylhs]);
  }

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
      yyerror (&yylloc, driver, YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyesa, &yyes, &yyes_capacity, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        if (yychar != YYEMPTY)
          YY_LAC_ESTABLISH;
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
        yyerror (&yylloc, driver, yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }

  yyerror_range[1] = yylloc;

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
                      yytoken, &yylval, &yylloc, driver);
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

      yyerror_range[1] = *yylsp;
      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp, yylsp, driver);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  /* If the stack popping above didn't lose the initial context for the
     current lookahead token, the shift below will for sure.  */
  YY_LAC_DISCARD ("error recovery");

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  yyerror_range[2] = yylloc;
  /* Using YYLLOC is tempting, but would change the location of
     the lookahead.  YYLOC is available though.  */
  YYLLOC_DEFAULT (yyloc, yyerror_range, 2);
  *++yylsp = yyloc;

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

#if 1
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (&yylloc, driver, YY_("memory exhausted"));
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
                  yytoken, &yylval, &yylloc, driver);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp, yylsp, driver);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  if (yyes != yyesa)
    YYSTACK_FREE (yyes);
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
