/* A Bison parser, made by GNU Bison 3.2.  */

/* Bison interface for Yacc-like parsers in C

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

/* Undocumented macros, especially those whose name start with YY_,
   are private implementation details.  Do not rely on them.  */

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
#line 37 "../../projects/generate_confs/parser.ypp" /* yacc.c:1906  */

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

#line 388 "generate_confs/parser.hpp" /* yacc.c:1906  */
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
