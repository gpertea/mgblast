/* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================*/

/** @file composition_adjustment.c
 * Highest level functions to solve the optimization problem for
 * compositional score matrix adjustment.
 *
 * @author Yi-Kuo Yu, Alejandro Schaffer, E. Michael Gertz
 */
#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] =
    "$Id: composition_adjustment.c,v 1.17 2006/04/28 15:07:47 gertz Exp $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include <limits.h>
#include <assert.h>
#include <algo/blast/core/ncbi_std.h>
#include <algo/blast/composition_adjustment/composition_constants.h>
#include <algo/blast/composition_adjustment/composition_adjustment.h>
#include <algo/blast/composition_adjustment/matrix_frequency_data.h>
#include <algo/blast/composition_adjustment/nlm_linear_algebra.h>
#include <algo/blast/composition_adjustment/optimize_target_freq.h>

/**positions of true characters in protein alphabet*/
static int trueCharPositions[COMPO_NUM_TRUE_AA] =
{1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22};

/**
 * conversion from 26 letter NCBIstdaa alphabet to 20 letter order
 * for true amino acids: ARNDCQEGHILKMFPSTWYV.  This order is
 * alphabetical in the standard three-letter abbreviation of each
 * amino acid */
static int alphaConvert[COMPO_PROTEIN_ALPHABET] =
  {(-1), 0, (-1),  4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15,
   16, 19, 17, (-1), 18, (-1), (-1), (-1)};


/**
 * Desired margin between an end of region used for computing a
 * composition, and the nearest StopChar; the desired margin may
 * not be attained. */
static const int kCompositionMargin = 20;

#define SCORE_BOUND            0.0000000001 /**< average scores below
                                                 -SCORE_BOUND are considered
                                                 effectively nonnegative, and
                                                 Newton's method will
                                                 will terminate */
#define LAMBDA_ITERATION_LIMIT 100          /**< iteration limit for Newton's
                                                 method. */
#define LAMBDA_ERROR_TOLERANCE 0.0000001    /**< bound on error for estimating
                                                 lambda */
#define LAMBDA_FUNCTION_TOLERANCE 0000001   /**< bound on the difference 
                                                 between the expected value
                                                 of exp(lambda x S) and 1 */

/** bound on error for Newton's method */
static const double kCompoAdjustErrTolerance = 0.00000001;
/** iteration limit for Newton's method */
static const int kCompoAdjustIterationLimit = 2000;
/** relative entropy of BLOSUM62 */
static const double kFixedReBlosum62 = 0.44;


/* Documented in composition_adjustment.h. */
void
Blast_ApplyPseudocounts(double * probs,
                    int number_of_observations,
                    const double * background_probs,
                    int pseudocounts)
{
    int i;                 /* loop index */
    double weight;         /* weight assigned to pseudocounts */
    double sum;            /* sum of the observed frequencies */
    /* pseudocounts as a double */
    double dpseudocounts = pseudocounts;

    /* Normalize probabilities */
    sum = 0.0;
    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        sum += probs[i];
    }
    if (sum == 0.0) {  /* Can't normalize a zero vector */
        sum = 1.0;
    }
    weight = dpseudocounts / (number_of_observations + dpseudocounts);
    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        probs[i] = (1.0 - weight) * probs[i] / sum
            + weight * background_probs[i];
    }
}


/* Documented in composition_adjustment.h. */
double
Blast_GetRelativeEntropy(const double A[], const double B[])
{
    int i;                 /* loop index over letters */
    double temp;           /* intermediate term */
    double value = 0.0;    /* square of relative entropy */

    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        temp = (A[i] + B[i]) / 2;
        if (temp > 0) {
            if (A[i] > 0) {
                value += A[i] * log(A[i] / temp) / 2;
            }
            if (B[i] > 0) {
                value += B[i] * log(B[i] / temp) / 2;
            }
        }
    }
    if (value < 0) {             /* must be numerical rounding error */
        value = 0;
    }
    return sqrt(value);
}



/* Blast_CalcLambdaFullPrecision -- interface documented in
 * composition_adjustment.h.
 *
 * If the average score for a composition is negative, and the maximum
 * score that occurs with nonzero probability is positive, then
 * statistical parameter Lambda exists and is the unique, positive
 * solution to
 *
 *    phi(lambda) = sum_{i,j} P_1(i) P_2(j) exp(S_{ij} lambda) - 1 = 0,
 *
 * where S_{ij} is the matrix "score" and P_1 and P_2 are row_probs and
 * col_probs respectively.
 *
 * It is simpler to solve this problem in x = exp(-lambda) than it is
 * to solve it in lambda, because we know that for x, a solution lies
 * in [0,1].  Furthermore, if M is the largest S_{ij} so that P_1(i)
 * and P_2(j) are nonzero, then the function
 *
 *    f(x) = -x^M - sum_{i,j} P_1(i) P_2(j) x^{M - S_{ij}},
 *
 * obtained by multiplying phi(lambda) by x^M, is well behaved in
 * (0,1] -- if the scores are integers, it is a polynomial.  Since x = 0
 * is not a solution, x solves f(x) = 0 in [0,1), if and only if
 * lambda = -ln(x) is a positive solution to phi(lambda).  Therefore,
 * we may define a safeguarded Newton iteration to find a solution of
 * f(x) = 0.
 *
 * For the most part, this is a standard safeguarded Newton iteration:
 * define an interval of uncertainty [a,b] with f(a) > 0 and f(b) < 0
 * (except for the initial value b = 1, where f(b) = 0); evaluate the
 * function and use the sign of that value to shrink the interval of
 * uncertainty; compute a Newton step; and if the Newton step suggests
 * a point outside the interval of uncertainty or fails to decrease
 * the function sufficiently, then bisect.  There are three further
 * details needed to understand the algorithm:
 *
 * 1)  If y the unique solution in [0,1], then f is positive to the left of
 *     y, and negative to the right.  Therefore, we may determine whether
 *     the Newton step -f(x)/f'(x) is moving toward, or away from, y by
 *     examining the sign of f'(x).  If f'(x) >= 0, we bisect, instead
 *     of taking the Newton step.
 * 2)  There is a neighborhood around x = 1 for which f'(x) >= 0,
 *     so (1) prevents convergence to x = 1.
 * 3)  Conditions like  fabs(p) < lambda_tolerance * x * (1-x) are used in
 *     convergence criteria because these values translate to a bound
 *     on the relative error in lambda.  This is proved in the
 *     "Blast Scoring Parameters" document that accompanies the BLAST
 *     code.
 *
 * We have observed that in typical cases the safeguarded Newton
 * iteration on f(x) requires half the iterations of a Newton
 * iteration on phi(lambda). More importantly, the iteration on f(x)
 * is robust and doesn't overflow; defining a robust safeguarded
 * Newton iteration on phi(lambda) that cannot converge to zero and
 * that is protected against overflow is more difficult.  So (despite
 * the length of this comment) the Newton iteration on f(x) is the
 * simpler solution.
 */
void
Blast_CalcLambdaFullPrecision(double * plambda, int *piterations,
                              double ** score, int alphsize,
                              const double row_prob[],
                              const double col_prob[],
                              double lambda_tolerance,
                              double function_tolerance,
                              int max_iterations)
{
    double f = 4;               /* The current function value; initially
                                   set to a value greater than any possible
                                   real value of f */
    double left = 0, right = 1; /* (left, right) is an interval containing
                                   a solution */
    double x = 0.367879441171;  /* The current iterate; initially exp(-1) */
    int is_newton = 0;          /* true if the last iteration was a Newton
                                   step; initially false */
    int i, j, k;                /* iteration indices */
    /* maximum score that occurs with nonzero probability */
    double max_score = COMPO_SCORE_MIN;
    /* average score */
    double avg_score = 0.0;

    /* Find the maximum score with nonzero probability */
    for (i = 0;  i < alphsize;  i++) {
        if (row_prob[i] == 0.0) {
            continue;
        }
        for (j = 0;  j < alphsize;  j++) {
            if (col_prob[j] == 0.0) {
                continue;
            }
            if (max_score < score[i][j]) {
                max_score = score[i][j];
            }
            avg_score += row_prob[i] * col_prob[j] * score[i][j];
        }
    }
    if (max_score <= 0.0 || avg_score >= 0) { 
        /* The iteration cannot converge if max_score is nonpositive
         * or the average score is nonnegative; lambda doesn't exist */
        *piterations = max_iterations;
        *plambda = -1.0;
        return;
    }
    for (k = 0;  k < max_iterations;  k++) {
        double slope;               /* slope of f at x */
        double fold = f;            /* previous value of f */
        double x_pow_max_score;     /* x raised to the power max_score */
        double lambda = -log(x);    /* the iterate lambda, see above */
        int was_newton = is_newton; /* true if the previous iteration
                                       was a Newton step; instead of a
                                       bisection step */
        /* Evaluate the function and its derivative */
        x_pow_max_score = exp(-max_score * lambda);
        f = -x_pow_max_score;
        slope = max_score * f / x;
        for (i = 0;  i < alphsize;  i++) {
            if (row_prob[i] == 0.0) {
                continue;
            }
            for (j = 0;  j < alphsize;  j++) {
                double ff;  /* a term in the sum used to compute f */

               if (col_prob[j] == 0.0) {
                   continue;
               }
               if (max_score != score[i][j]) {
                   double diff_score = max_score - score[i][j];

                   ff = row_prob[i] * col_prob[j] * exp(-lambda * diff_score);
                   slope += diff_score * ff / x;
               } else {
                   ff = row_prob[i] * col_prob[j];
               }
               f += ff;
            }
        }
        /* Finished evaluating the function and its derivative */
        if (f > 0) {
            left = x; /* move the left endpoint */
        } else if (f < 0) {
            right = x; /* move the right endpoint */
        } else { /* f == 0 */
            break; /* x is an exact solution */
        }
        if (right - left <= 2 * left * (1 - right) * lambda_tolerance &&
            fabs(f/x_pow_max_score) <= function_tolerance) {
            /* The midpoint of the interval converged */
            x = (left + right) / 2;
            break;
        }
        if ((was_newton && fabs(f) > .9 * fabs(fold))
            /* if the previous iteration was a Newton step but didn't
             * decrease f sufficiently; or */
             || slope >= 0
             /* if a Newton step will move us away from the desired solution */
        ) {/* then */
            x = (left + right)/2;  /* bisect */
        } else {
            double p = -f/slope;   /* The Newton step */
            double y = x + p;      /* The proposed next iterate */
            if (y <= left || y >= right) { /* The proposed iterate is
                                              not in (left,right) */
                x = (left + right)/2;  /* bisect */
            } else {/* The proposed iterate is in (left,right). Accept it. */
                is_newton = 1;
                x = y;
                if (fabs(p) <= lambda_tolerance * x * (1-x) &&
                    fabs(f/x_pow_max_score) <= function_tolerance) break;
            }
        }
    }  /* End for all iterations k */
    *plambda = -log(x);
    *piterations = k;
}


/* Documented in composition_adjustment.h. */
double
Blast_MatrixEntropy(double ** matrix, int alphsize, const double row_prob[],
                    const double col_prob[], double Lambda)
{
    int i, j;
    double entropy = 0.0;
    for (i = 0;  i < alphsize;  i++) {
        for (j = 0;  j < alphsize;  j++) {
            /* the score at (i,j), rescaled to nats */
            double nat_score = Lambda * matrix[i][j];
            entropy += nat_score * exp(nat_score) * row_prob[i] * col_prob[j];
        }
    }
    return entropy;
}



/* Documented in composition_adjustment.h. */
double
Blast_TargetFreqEntropy(double ** target_freq)
{
    int i, j;          /* iteration indices */
    double entropy;    /* the entropy to be returned */
    /* Row probabilities consistent with the target frequencies */
    double row_prob[COMPO_NUM_TRUE_AA] = {0,};
    double col_prob[COMPO_NUM_TRUE_AA] = {0,};

    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        for (j = 0;  j < COMPO_NUM_TRUE_AA;  j++) {
            row_prob[i] += target_freq[i][j];
            col_prob[j] += target_freq[i][j];
        }
    }
    entropy = 0.0;
    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        for (j = 0;  j < COMPO_NUM_TRUE_AA;  j++) {
            double freq = target_freq[i][j];
            entropy += freq * log(freq / row_prob[i] / col_prob[j]);
        }
    }
    return entropy;
}


/* Documented in composition_adjustment.h. */
void
Blast_CalcFreqRatios(double ** ratios, int alphsize,
                     double row_prob[], double col_prob[])
{
    int i, j;
    for (i = 0;  i < alphsize;  i++) {
        if (row_prob[i] > 0) {
            for (j = 0;  j < alphsize;  j++) {
                if (col_prob[j] > 0) {
                    ratios[i][j] /= (row_prob[i] * col_prob[j]);
                }
            }
        }
    }
}


/* Documented in composition_adjustment.h. */
void
Blast_FreqRatioToScore(double ** matrix, int rows, int cols, double Lambda)
{
    int i;
    for (i = 0;  i < rows;  i++) {
        int j;
        for (j = 0;  j < cols;  j++) {
            if (0.0 == matrix[i][j]) {
                matrix[i][j] = COMPO_SCORE_MIN;
            } else {
                matrix[i][j] = log(matrix[i][j])/Lambda;
            }
        }
    }
}


/**
 * Convert letter probabilities from the NCBIstdaa alphabet to
 * a 20 letter ARND... amino acid alphabet. (@see alphaConvert)
 *
 * @param outputLetterProbs   the ARND letter probabilities [out]
 * @param inputLetterProbs    the NCBIstdaa letter probabilities [in]
 */
static void
s_GatherLetterProbs(double * outputLetterProbs,
                    const double * inputLetterProbs)
{
    int c; /*index over characters*/

    for (c = 0;  c < COMPO_PROTEIN_ALPHABET;  c++) {
        if ((-1) != alphaConvert[c]) {
            outputLetterProbs[alphaConvert[c]] = inputLetterProbs[c];
        }
    }
}

/**
 * Convert letter probabilities from a a 20 letter ARND... amino acid
 * alphabet to the NCBIstdaa alphabet, (@see alphaConvert) setting
 * probabilities for nonstandard characters to zero.
 *
 * @param std_probs   the NCBIstdaa letter probabilities [out]
 * @param probs       the ARND letter probabilities [in]
 */
static void
s_UnpackLetterProbs(double std_probs[], const double probs[])
{
    int c; /*index over characters*/

    for (c = 0;  c < COMPO_PROTEIN_ALPHABET;  c++) {
        if ((-1) != alphaConvert[c]) {
            std_probs[c] = probs[alphaConvert[c]];
        } else {
            std_probs[c] = 0.0;
        }
    }
}


/**
 * Set the probabilities for the two-character ambiguity characters.
 * @param probs     the probabilities for the NCBIstdaa alphabet.
 *                  On entry, values for the standard amino acids must
 *                  be set; on exit, values for the two-character
 *                  ambiguities will also be set.
 */
static void
s_SetPairAmbigProbsToSum(double probs[])
{
    probs[eBchar] = probs[eDchar] + probs[eNchar];
    probs[eZchar] = probs[eEchar] + probs[eQchar];
}


/* Documented in composition_adjustment.h. */
int
Blast_EntropyOldFreqNewContext(double * entropy,
                               double * Lambda,
                               int * iter_count,
                               double ** target_freq,
                               const double row_prob[],
                               const double col_prob[])
{
    /* iteration indices */
    int i, j;
    /* Status flag; will be set to zero on success */
    int status = 1;
    /* A matrix of scores in the context constistent with the target 
     * frequencies */
    double ** scores;
    /* Row and column probabilities consistent with the target
     * frequencies; the old context */
    double old_col_prob[COMPO_NUM_TRUE_AA] = {0.0,};
    double old_row_prob[COMPO_NUM_TRUE_AA] = {0.0,};

    *entropy = 0;
    status = 1;

    /* Calculate the matrix "scores" from the target frequencies */
    scores = Nlm_DenseMatrixNew(COMPO_NUM_TRUE_AA, COMPO_NUM_TRUE_AA);
    if (scores == NULL) {
        return -1;
    }
    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        for (j = 0;  j < COMPO_NUM_TRUE_AA; j++) {
            old_row_prob[i] += target_freq[i][j];
            old_col_prob[j] += target_freq[i][j];
        }
    }
    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        memcpy(scores[i], target_freq[i], COMPO_NUM_TRUE_AA * sizeof(double));
    }
    Blast_CalcFreqRatios(scores, COMPO_NUM_TRUE_AA,
                         old_row_prob, old_col_prob);
    Blast_FreqRatioToScore(scores, COMPO_NUM_TRUE_AA, COMPO_NUM_TRUE_AA, 1.0);
    /* Finished calculating the matrix "scores" */

    Blast_CalcLambdaFullPrecision(Lambda, iter_count, scores,
                                  COMPO_NUM_TRUE_AA, row_prob,
                                  col_prob, LAMBDA_ERROR_TOLERANCE,
                                  LAMBDA_FUNCTION_TOLERANCE,
                                  LAMBDA_ITERATION_LIMIT);
    if (*iter_count <  LAMBDA_ITERATION_LIMIT) {
        *entropy = Blast_MatrixEntropy(scores, COMPO_NUM_TRUE_AA,
                                       row_prob, col_prob, *Lambda);
        status = 0;
    }
    Nlm_DenseMatrixFree(&scores);
    return status;
}


/**
 * Convert a matrix of target frequencies for the ARND alphabet of
 * true amino acids to a set of target frequencies for the NCBIstdaa
 * alphabet, filling in value for the two-character ambiguities (but
 * not X).
 *
 * @param StdFreq      frequencies in the NCBIstdaa alphabet [output]
 * @param freq         frequencies in the ARND alphabet [input]
 */
static void
s_TrueAaToStdTargetFreqs(double ** StdFreq, double ** freq)
{
    /* Note I'm using a rough convention for this routine that uppercase
     * letters refer to quantities in the standard (larger) alphabet
     * and lowercase letters refer to the true amino acid (smaller)
     * alphabet.
     */
    /* Shorter names for the sizes of the two alphabets */
    const int small_alphsize = COMPO_NUM_TRUE_AA;
    const int StdAlphsize = COMPO_PROTEIN_ALPHABET;
    int A, B;          /* characters in the std (big) alphabet */
    int a, b;          /* characters in the small alphabet */
    double sum;        /* sum of values in target_freq; used to normalize */
    sum = 0.0;
    for (a = 0;  a < small_alphsize;  a++) {
        for (b = 0;  b < small_alphsize;  b++) {
            sum +=  freq[a][b];
        }
    }
    for (A = 0;  A < StdAlphsize;  A++) {
        /* for all rows */
        if (alphaConvert[A] < 0) {
            /* the row corresponds to a nonstandard reside */
            for (B = 0;  B < StdAlphsize;  B++) {
                StdFreq[A][B] = 0.0;
            }
        } else {
            /* the row corresponds to a standard reside */
            a = alphaConvert[A];

            for (B = 0;  B < StdAlphsize;  B++) {
                /* for all columns */
                if (alphaConvert[B] < 0) {
                    /* the column corresponds to a nonstandard reside */
                    StdFreq[A][B] = 0.0;
                } else {
                    /* the column corresponds to a standard reside */
                    b = alphaConvert[B];
                    StdFreq[A][B] = freq[a][b] / sum;
                }
            }
            /* Set values for two-character ambiguities */
            StdFreq[A][eBchar] = StdFreq[A][eDchar] + StdFreq[A][eNchar];
            StdFreq[A][eZchar] = StdFreq[A][eEchar] + StdFreq[A][eQchar];
        }
    }
    /* Add rows to set values for two-character ambiguities */
    memcpy(StdFreq[eBchar], StdFreq[eDchar], StdAlphsize * sizeof(double));
    Nlm_AddVectors(StdFreq[eBchar], StdAlphsize, 1.0, StdFreq[eNchar]);

    memcpy(StdFreq[eZchar], StdFreq[eEchar], StdAlphsize * sizeof(double));
    Nlm_AddVectors(StdFreq[eZchar], StdAlphsize, 1.0, StdFreq[eQchar]);
}


/**
 * Calculate values for the X ambiguity score, give an vector of scores and
 * a set of character probabilities.
 * @param M            a vector of scores
 * @param incM         stride between elements of M; used to pass a column
 *                     of a score matrix to this routine by setting incM
 *                     to the number of elements in a column.
 * @param probs        background probabilities for a sequence.
 */
static double
s_CalcXScore(double * M, int incM, const double probs[])
{
    int j;                   /* iteration index */
    double score_iX = 0.0;   /* score of character i substituted by X */

    for (j = 0;  j < COMPO_PROTEIN_ALPHABET;  j++) {
        if (alphaConvert[j] >= 0) {
            /* If the column corresponds to a true amino acid */
            score_iX += M[j * incM] * probs[j];
        }
    }
    return score_iX <= -1.0 ? score_iX : -1.0;
}


/**
 * Given a standard matrix with the scores for true amino acid, and
 * pairwise ambiguity scores set, calculate and set the scores for the
 * X and U (Selenocystiene) characters.
 *
 * @param M             a scoring matrix
 * @param row_probs     character frequencies for the sequence corresponding
 *                      to the rows of M.
 * @param col_probs     character frequencies for the sequence corresponding
 *                      to the columns of M.
 */
static void
s_SetXUScores(double ** M, const double row_probs[], const double col_probs[])
{
    int i;                      /* iteration index */
    double score_XX = 0.0;      /* score of matching an X to an X */
    /* A shorter name for COMPO_PROTEIN_ALPHABET */
    const int alphsize = COMPO_PROTEIN_ALPHABET;

    for (i = 0;  i < alphsize;  i++) {
        if (alphaConvert[i] >= 0) {
            M[i][eXchar] = s_CalcXScore(M[i], 1, col_probs);
            score_XX += M[i][eXchar] * row_probs[i];

            M[eXchar][i] = s_CalcXScore(&M[0][i], alphsize, row_probs);
        }
    }
    M[eXchar][eXchar] = score_XX;

    /* Set X scores for pairwise ambiguity characters */
    M[eBchar][eXchar] = s_CalcXScore(M[eBchar], 1, col_probs);
    M[eXchar][eBchar] = s_CalcXScore(&M[0][eBchar], alphsize, row_probs);

    M[eZchar][eXchar] = s_CalcXScore(M[eZchar], 1, col_probs);
    M[eXchar][eZchar] = s_CalcXScore(&M[0][eZchar], alphsize, row_probs);
    /* Copy X scores to U */
    memcpy(M[eSelenocysteine], M[eXchar], alphsize * sizeof(double));
    for (i = 0;  i < alphsize;  i++) {
        M[i][eSelenocysteine] = M[i][eXchar];
    }
}


/** Return the nearest integer to x. */
static long Nint(double x)
{
    x += (x >= 0. ? 0.5 : -0.5);
    return (long)x;
}


/**
 * Round a floating point matrix of scores to produce an integer matrix of
 * scores.
 *
 * @param matrix             the matrix of integer valued scores [out]
 * @param floatScoreMatrix   the matrix of floating point valued
 *                           scores [in]
 * @param numPositions       the number of rows of the matrices.
 */
static void
s_RoundScoreMatrix(int **matrix, int numPositions,
                   double **floatScoreMatrix)
{
    int p, c; /*indices over positions and characters*/

    for (p = 0;  p < numPositions;  p++) {
        for (c = 0;  c < COMPO_PROTEIN_ALPHABET;  c++) {
            if (floatScoreMatrix[p][c] < INT_MIN) {
                matrix[p][c] = INT_MIN;
            } else {
                matrix[p][c] = Nint(floatScoreMatrix[p][c]);
            }
        }
    }
}


/**
 * Given a set of target frequencies and two sets of character
 * probabilities for the true amino acids in the ARND alphabet,
 * calculate a scoring matrix that has valid entries for all
 * characters in the NCBIstdaa amino acid alphabet.
 *
 * @param Matrix        the newly computed matrix
 * @param target_freq   target frequencies for true amino acids (20x20)
 * @param StartMatrix   a matrix containing values for the stop character
 * @param row_prob      probabilities of true amino acids in the sequence
 *                      corresponding to the rows of matrix (length = 20)
 * @param col_prob      probabilities of true amino acids in the sequence
 *                      corresponding to the columns of matrix (length = 20)
 * @param Lambda        the desired scale of this matrix
 */
static int
s_ScoresStdAlphabet(int ** Matrix, double ** target_freq, int ** StartMatrix,
                    const double row_prob[], const double col_prob[],
                    double Lambda)
{
    /* Note: I'm using a rough convention for this routine that uppercase
     * letters refer to quantities in the standard (larger) alphabet
     * and lowercase letters refer to the true amino acid (smaller)
     * alphabet.
     */
    int i;
    /* row and column probabilities in the NCBIstdaa alphabet */
    double RowProb[COMPO_PROTEIN_ALPHABET];
    double ColProb[COMPO_PROTEIN_ALPHABET];
    /* A double precision score matrix */
    double ** Scores = Nlm_DenseMatrixNew(COMPO_PROTEIN_ALPHABET,
                                          COMPO_PROTEIN_ALPHABET);
    if (Scores == NULL) {
        return -1;
    }
    s_UnpackLetterProbs(RowProb, row_prob);
    s_SetPairAmbigProbsToSum(RowProb);

    s_UnpackLetterProbs(ColProb, col_prob);
    s_SetPairAmbigProbsToSum(ColProb);

    s_TrueAaToStdTargetFreqs(Scores, target_freq);
    Blast_CalcFreqRatios(Scores, COMPO_PROTEIN_ALPHABET, RowProb, ColProb);
    Blast_FreqRatioToScore(Scores, COMPO_PROTEIN_ALPHABET,
                           COMPO_PROTEIN_ALPHABET, Lambda);
    s_SetXUScores(Scores, RowProb, ColProb);

    s_RoundScoreMatrix(Matrix, COMPO_PROTEIN_ALPHABET, Scores);
    Nlm_DenseMatrixFree(&Scores);

    for (i = 0;  i < COMPO_PROTEIN_ALPHABET;  i++) {
        Matrix[i][eStopChar] = StartMatrix[i][eStopChar];
        Matrix[eStopChar][i] = StartMatrix[eStopChar][i];
    }
    return 0;
}


/**
 * Find the range of scores contained in an scoring matrix.
 * @param obs_min    smallest value in the matrix
 * @param obs_max    largest value in the matrix
 * @param matrix     a matrix with COMPO_NUM_TRUE_AA columns
 * @param rows       number of rows in the matrix
 */
static void s_GetScoreRange(int * obs_min, int * obs_max,
                            int ** matrix, int rows)
{
    int aa;                    /* index of an amino-acid in the 20
                                  letter alphabet */
    int irow, jcol;            /* matrix row and column indices */
    int minScore, maxScore;    /* largest and smallest observed scores */

    minScore = maxScore = 0;
    for (irow = 0;  irow < rows;  irow++) {
        for (aa = 0;  aa < COMPO_NUM_TRUE_AA;  aa++) {
            jcol = trueCharPositions[aa];
            if (matrix[irow][jcol] < minScore &&
                matrix[irow][jcol] > COMPO_SCORE_MIN)
                minScore = matrix[irow][jcol];
            if (matrix[irow][jcol] > maxScore)
                maxScore = matrix[irow][jcol];
        }
    }
    *obs_min = minScore;
    *obs_max = maxScore;
}


/**
 * Compute the score probabilities for a given amino acid substitution matrix
 * in the context of given query and subject amino acid frequencies.
 *
 * @param *obs_min          the smallest score in the score matrix [out]
 * @param *obs_max          the largest score in the score matrix [out]
 * @param *scoreProb        the new array, of length (*obs_max - *obs_min + 1),
 *                          of score probabilities, where (*scoreProb)[0] is
 *                          the probability for score *obs_min.
 * @param matrix            a amino-acid substitution matrix (not
 *                          position-specific)
 * @param subjectProbArray  is an array containing the probability of
 *                          occurrence of each residue in the subject
 * @param queryProbArray    is an array containing the probability of
 *                          occurrence of each residue in the query
 * @param scoreProb         is an array of probabilities for each score
 *                          that is to be used as a field in return_sfp
 * @return 0 on success, -1 on out-of-memory
 */
static int
s_GetMatrixScoreProbs(double **scoreProb, int * obs_min, int * obs_max,
                      int **matrix, const double *subjectProbArray,
                      const double *queryProbArray)
{
    int aa;          /* index of an amino-acid in the 20 letter
                        alphabet */
    int irow, jcol;  /* matrix row and column indices */
    double * sprob;  /* a pointer to the element of the score
                        probabilities array that represents the
                        probability of the score 0*/
    int minScore;    /* smallest score in matrix; the same value as
                        (*obs_min). */
    int range;       /* the range of scores in the matrix */

    s_GetScoreRange(obs_min, obs_max, matrix, COMPO_PROTEIN_ALPHABET);
    minScore = *obs_min;
    range = *obs_max - *obs_min + 1;
    *scoreProb = calloc(range, sizeof(double));
    if (*scoreProb == NULL) {
        return -1;
    }
    sprob = &((*scoreProb)[-(*obs_min)]); /*center around 0*/
    for (irow = 0;  irow < COMPO_PROTEIN_ALPHABET;  irow++) {
        for (aa = 0;  aa < COMPO_NUM_TRUE_AA;  aa++) {
            jcol = trueCharPositions[aa];
            if (matrix[irow][jcol] >= minScore) {
                sprob[matrix[irow][jcol]] +=
                    (queryProbArray[irow] * subjectProbArray[jcol]);
            }
        }
    }
    return 0;
}


/**
 * Compute the score probabilities for a given amino acid position-specific
 * substitution matrix in the context of a given set of subject amino
 * acid frequencies.
 *
 * @param *obs_min          the smallest score in the score matrix [out]
 * @param *obs_max          the largest score in the score matrix [out]
 * @param *scoreProb        the new array, of length (*obs_max - *obs_min + 1),
 *                          of score probabilities, where (*scoreProb)[0] is
 *                          the probability for score *obs_min.
 * @param matrix            a position-specific amino-acid substitution matrix.
 * @param rows              the number of rows in matrix.
 * @param subjectProbArray  is an array containing the probability of
 *                          occurrence of each residue in the subject
 * @return 0 on success, -1 on out-of-memory
 */
static int
s_GetPssmScoreProbs(double ** scoreProb, int * obs_min, int * obs_max,
                    int **matrix, int rows,
                    const double *subjectProbArray)
{
    int aa;            /* index of an amino-acid in the 20 letter
                          alphabet */
    int irow, jcol;    /* matrix row and column indices */
    double onePosFrac; /* matrix length as a double*/
    double * sprob;    /* pointer to the element of the score
                        * probabilities array the represents the
                        * probability of zero */
    int minScore;      /* smallest score in matrix; the same value as
                          (*obs_min). */
    int range;         /* the range of scores in the matrix */

    s_GetScoreRange(obs_min, obs_max, matrix, rows);
    minScore = *obs_min;
    range = *obs_max - *obs_min + 1;
    *scoreProb = calloc(range, sizeof(double));
    if (*scoreProb == NULL) {
        return -1;
    }
    sprob = &((*scoreProb)[-(*obs_min)]); /*center around 0*/
    onePosFrac = 1.0/ ((double) rows);
    for (irow = 0;  irow < rows;  irow++) {
        for (aa = 0;  aa < COMPO_NUM_TRUE_AA;  aa++) {
            jcol = trueCharPositions[aa];
            if (matrix[irow][jcol] >= minScore) {
                sprob[matrix[irow][jcol]] +=
                    onePosFrac * subjectProbArray[jcol];
            }
        }
    }
    return 0;
}


/* Documented in composition_adjustment.h. */
void
Blast_Int4MatrixFromFreq(Int4 **matrix, int alphsize,
                         double ** freq, double Lambda)
{
    /* TODO: Eliminate this routine, or change its API.  The void return
     * value means that the routine cannot fail (and so cannot allocate
     * memory) but alphsize suggests that it can use an arbitrarily large
     * alphabet, while in truth the largest alphabet it can handle is
     * COMPO_PROTEIN_ALPHABET.
     *
     * The routine exists in its current form because API changes
     * to this library are very disruptive.  It should be changed the next
     * time an API change is inevitable */

    double dMatrixStore[COMPO_PROTEIN_ALPHABET];
    double * dMatrix[1];
    int i;
    
    dMatrix[0] = dMatrixStore;

    assert(alphsize <= COMPO_PROTEIN_ALPHABET);

    for (i = 0;  i < alphsize;  i++) {
        memcpy(dMatrix[0], freq[i], alphsize * sizeof(double));
        Blast_FreqRatioToScore(dMatrix, 1, COMPO_PROTEIN_ALPHABET, Lambda);
        s_RoundScoreMatrix(&matrix[i], 1, dMatrix);
    }
}


/* Documented in composition_adjustment.h. */
void Blast_MatrixInfoFree(Blast_MatrixInfo ** ss)
{
    if (*ss != NULL) {
        free((*ss)->matrixName);
        Nlm_Int4MatrixFree(&(*ss)->startMatrix);
        Nlm_DenseMatrixFree(&(*ss)->startFreqRatios);
        free(*ss);
        *ss = NULL;
    }
}


/* Documented in composition_adjustment.h. */
Blast_MatrixInfo *
Blast_MatrixInfoNew(int rows, int positionBased)
{
    int i;       /* loop index */
    Blast_MatrixInfo * ss = malloc(sizeof(Blast_MatrixInfo));
    if (ss != NULL) {
        ss->rows = rows;
        ss->positionBased = positionBased;

        ss->matrixName = NULL;
        ss->startMatrix = NULL;
        ss->startFreqRatios = NULL;

        ss->startMatrix  = Nlm_Int4MatrixNew(rows + 1, COMPO_PROTEIN_ALPHABET);
        if (ss->startMatrix == NULL)
            goto error_return;
        ss->startFreqRatios = Nlm_DenseMatrixNew(rows + 1,
                                                 COMPO_PROTEIN_ALPHABET);
        if (ss->startFreqRatios == NULL)
            goto error_return;
        for (i = 0;  i < COMPO_PROTEIN_ALPHABET;  i++) {
            ss->startMatrix[rows][i] = COMPO_SCORE_MIN;
            ss->startFreqRatios[rows][i] = (double) COMPO_SCORE_MIN;
        }

    }
    goto normal_return;
error_return:
    Blast_MatrixInfoFree(&ss);
normal_return:
    return ss;
}


/**
 * Fill in all scores for a PSSM at a given scale Lambda.
 *
 * @param matrix       the newly computed matrix [output]
 * @param rows         number of positions (rows) in the PSSM
 * @param freq_ratios  frequency ratios defining the PSSM
 * @param start_matrix an existing PSSM; used to set values for the
 *                     stop character.
 * @param col_prob     letter probabilities
 * @param Lambda       scale of the new matrix
 */
static void
s_ScalePSSM(int **matrix, int rows, double ** freq_ratios,
            int ** start_matrix, const double col_prob[], double Lambda)
{
    int p;          /* index over matrix rows */
    /* A row of scores corresponding to one position in the PSSM */
    double row[COMPO_PROTEIN_ALPHABET];
    /* A matrix with one row */
    double * row_matrix[1];

    row_matrix[0] = row;

    for (p = 0;  p < rows;  p++) {
        memcpy(row, freq_ratios[p], sizeof(row));

        Blast_FreqRatioToScore(row_matrix, 1, COMPO_PROTEIN_ALPHABET, Lambda);
        row[eXchar] = s_CalcXScore(row, 1, col_prob);
        row[eSelenocysteine] = row[eXchar];
        s_RoundScoreMatrix(&matrix[p], 1, row_matrix);

        matrix[p][eStopChar] = start_matrix[p][eStopChar];
    }
}


/**
 * Fill in all scores for a scoring matrix for the NCBIstdaa alphabet
 * at a given scale Lambda.
 *
 * @param matrix       the newly computed matrix [output]
 * @param freq_ratios  frequency ratios defining the PSSM
 * @param start_matrix an existing matrix; used to set values for the
 *                     stop character.
 * @param row_prob     letter probabilities for the sequence
 *                     corresponding to the rows of the matrix.
 * @param col_prob     letter probabilities for the sequence
 *                     corresponding to the columns of the matrix.
 * @param Lambda       scale of the new matrix
 */
static int
s_ScaleSquareMatrix(int **matrix, double ** freq_ratios, int ** start_matrix,
                    const double row_prob[], const double col_prob[],
                    double Lambda)
{
    /* A shorter name for COMPO_PROTEIN_ALPHABET; the size of the alphabet */
    const int alphsize = COMPO_PROTEIN_ALPHABET;
    double ** scores;     /* a double precision matrix of scores */
    int i;                /* iteration index */

    scores = Nlm_DenseMatrixNew(alphsize, alphsize);
    if (scores == 0) return -1;

    for (i = 0;  i < alphsize;  i++) {
        memcpy(scores[i], freq_ratios[i], alphsize * sizeof(double));
    }
    Blast_FreqRatioToScore(scores, alphsize, alphsize, Lambda);
    s_SetXUScores(scores, row_prob, col_prob);
    s_RoundScoreMatrix(matrix, alphsize, scores);
    for (i = 0;  i < alphsize;  i++) {
        matrix[i][eStopChar] = start_matrix[i][eStopChar];
        matrix[eStopChar][i] = start_matrix[eStopChar][i];
    }
    Nlm_DenseMatrixFree(&scores);

    return 0;
}


/** LambdaRatioLowerBound is used when the expected score is too large
 * causing impalaKarlinLambdaNR to give a Lambda estimate that
 * is too small, or to fail entirely returning -1*/
#define LambdaRatioLowerBound 0.5


/* Documented in composition_adjustment.h. */
int
Blast_CompositionBasedStats(int ** matrix, double * LambdaRatio,
                            const Blast_MatrixInfo * ss,
                            const double queryProb[], const double resProb[],
                            double (*calc_lambda)(double*,int,int,double))
{
    double correctUngappedLambda; /* new value of ungapped lambda */
    int obs_min, obs_max;         /* smallest and largest score in the
                                     unscaled matrix */
    double *scoreArray;           /* an array of score probabilities */
    int out_of_memory;            /* status flag to indicate out of memory */

    if (ss->positionBased) {
        out_of_memory =
            s_GetPssmScoreProbs(&scoreArray, &obs_min, &obs_max,
                                ss->startMatrix, ss->rows, resProb);
    } else {
        out_of_memory =
            s_GetMatrixScoreProbs(&scoreArray, &obs_min, &obs_max,
                                  ss->startMatrix, resProb, queryProb);
    }
    if (out_of_memory)
        return -1;
    correctUngappedLambda =
        calc_lambda(scoreArray, obs_min, obs_max, ss->ungappedLambda);

    /* calc_lambda will return -1 in the case where the
     * expected score is >=0; however, because of the MAX statement 3
     * lines below, LambdaRatio should always be > 0; the succeeding
     * test is retained as a vestige, in case one wishes to remove the
     * MAX statement and allow LambdaRatio to take on the error value
     * -1 */
    *LambdaRatio = correctUngappedLambda / ss->ungappedLambda;
    *LambdaRatio = MIN(1, *LambdaRatio);
    *LambdaRatio = MAX(*LambdaRatio, LambdaRatioLowerBound);

    if (*LambdaRatio > 0) {
        double scaledLambda = ss->ungappedLambda/(*LambdaRatio);
        if (ss->positionBased) {
            s_ScalePSSM(matrix, ss->rows, ss->startFreqRatios,
                        ss->startMatrix, resProb, scaledLambda);
        } else {
            s_ScaleSquareMatrix(matrix, ss->startFreqRatios, ss->startMatrix,
                                queryProb, resProb, scaledLambda);
        }
    }
    free(scoreArray);

    return 0;
}


/* Documented in composition_adjustment.h. */
void
Blast_ReadAaComposition(Blast_AminoAcidComposition * composition,
                        const Uint1 * sequence, int length)
{
    int i; /* iteration index */

    /* fields of composition as local variables */
    int numTrueAminoAcids = 0;
    double * prob = composition->prob;

    for (i = 0;  i < COMPO_PROTEIN_ALPHABET;  i++) {
        prob[i] = 0.0;
    }
    for (i = 0;  i < length;  i++) {
        if (alphaConvert[sequence[i]] >= 0) {
            prob[sequence[i]]++;
            numTrueAminoAcids++;
        }
    }
    composition->numTrueAminoAcids = numTrueAminoAcids;
    if (numTrueAminoAcids > 0) {
        for (i = 0;  i < COMPO_PROTEIN_ALPHABET;  i++) {
            prob[i] /= numTrueAminoAcids;
        }
    }
}


/* Documented in composition_adjustment.h. */
void
Blast_GetCompositionRange(int * pleft, int * pright,
                          const Uint1 * subject_data, int length,
                          int start, int finish)
{
    int i;                /* iteration index */
    int left, right;

    left = start;
    /* Search leftward for a StopChar */
    for (i = left;  i > 0;  i--) {
        if (subject_data[i - 1] == eStopChar) {
            /* We have found a StopChar. Unless the StopChar is
             * too close to the start of the subject region of the
             * HSP, */
            if (i + kCompositionMargin < left) {
                /* reset the left endpoint. */
                left = i + kCompositionMargin;
            }
            break;
        }
    }
    if (i == 0) {
        /* No stop codon was found to the left. */
        left = 0;
    }
    right = finish;
    /* Search rightward for a StopChar */
    for (i = right;  i < length;  i++) {
        if (subject_data[i] == eStopChar) {
            /* We have found a StopChar. Unless the StopChar is
             * too close to the end of the subject region of the
             * HSP, */
            if (i - kCompositionMargin > right) {
                /* reset the right endpoint */
                right = i - kCompositionMargin;
            }
            break;
        }
    }
    if (i == length) {
        /* No stop codon was found to the right. */
        right = length;
    }
    *pleft = left; *pright = right;
}


/* Documented in composition_adjustment.h. */
void
Blast_CompositionWorkspaceFree(Blast_CompositionWorkspace ** pNRrecord)
{
    Blast_CompositionWorkspace * NRrecord = *pNRrecord;

    if (NRrecord != NULL) {
        free(NRrecord->first_standard_freq);
        free(NRrecord->second_standard_freq);

        Nlm_DenseMatrixFree(&NRrecord->mat_final);
        Nlm_DenseMatrixFree(&NRrecord->mat_b);

        free(NRrecord);
    }
    pNRrecord = NULL;
}


/* Documented in composition_adjustment.h. */
Blast_CompositionWorkspace * Blast_CompositionWorkspaceNew()
{
    Blast_CompositionWorkspace * NRrecord;        /* record to allocate
                                                    and return */
    int i;                     /* loop index */

    NRrecord = (Blast_CompositionWorkspace *)
        malloc(sizeof(Blast_CompositionWorkspace));
    if (NRrecord == NULL) goto error_return;

    NRrecord->first_standard_freq      = NULL;
    NRrecord->second_standard_freq     = NULL;
    NRrecord->mat_final                = NULL;
    NRrecord->mat_b                    = NULL;

    NRrecord->first_standard_freq =
        (double *) malloc(COMPO_NUM_TRUE_AA * sizeof(double));
    if (NRrecord->first_standard_freq == NULL) goto error_return;

    NRrecord->second_standard_freq =
        (double *) malloc(COMPO_NUM_TRUE_AA * sizeof(double));
    if (NRrecord->second_standard_freq == NULL) goto error_return;

    NRrecord->mat_final   = Nlm_DenseMatrixNew(COMPO_NUM_TRUE_AA,
                                               COMPO_NUM_TRUE_AA);
    if (NRrecord->mat_final == NULL) goto error_return;

    NRrecord->mat_b       = Nlm_DenseMatrixNew(COMPO_NUM_TRUE_AA,
                                               COMPO_NUM_TRUE_AA);
    if (NRrecord->mat_b == NULL) goto error_return;

    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        NRrecord->first_standard_freq[i] =
            NRrecord->second_standard_freq[i] = 0.0;
    }

    goto normal_return;
error_return:
    Blast_CompositionWorkspaceFree(&NRrecord);
normal_return:
    return NRrecord;
}


/* Documented in composition_adjustment.h. */
int
Blast_CompositionWorkspaceInit(Blast_CompositionWorkspace * NRrecord,
                               const char *matrixName)
{
    if (0 == Blast_GetJointProbsForMatrix(NRrecord->mat_b,
                                          NRrecord->first_standard_freq,
                                          NRrecord->second_standard_freq,
                                          matrixName)) {
        return 0;
    } else {
        fprintf(stderr,
                "Matrix %s not currently supported for RE based adjustment\n",
                matrixName);
        return -1;
    }
}


/* Documented in composition_adjustment.h. */
int
Blast_CompositionMatrixAdj(int ** matrix,
                           EMatrixAdjustRule matrix_adjust_rule,
                           int length1,
                           int length2,
                           const double * stdaa_row_probs,
                           const double * stdaa_col_probs,
                           int pseudocounts,
                           double specifiedRE,
                           Blast_CompositionWorkspace * NRrecord,
                           const Blast_MatrixInfo * matrixInfo)
{
    int i;                         /* loop indices */
    static int total_iterations = 0;   /* total iterations among all
                                          calls to
                                          compute_new_score_matrix */
    int new_iterations = 0;        /* number of iterations in the most
                                      recent call to
                                      compute_new_score_matrix */
    static int max_iterations = 0; /* maximum number of iterations
                                      observed in a call to
                                      compute_new_score_matrix */
    int status;                    /* status code for operations that may
                                      fail */
    double row_probs[COMPO_NUM_TRUE_AA];
    double col_probs[COMPO_NUM_TRUE_AA];
    double RE_final;

    /*Is the relative entropy constrained? Behaves as boolean for now*/
    int constrain_rel_entropy =
        eUnconstrainedRelEntropy != matrix_adjust_rule;

    s_GatherLetterProbs(row_probs, stdaa_row_probs);
    s_GatherLetterProbs(col_probs, stdaa_col_probs);

    switch (matrix_adjust_rule) {
    case eUnconstrainedRelEntropy:
        /* Initialize to a arbitrary value; it won't be used */
        RE_final = 0.0;
        break;
    case eRelEntropyOldMatrixNewContext:
        {
            double entropy, Lambda;
            int iter_count;
            status = Blast_EntropyOldFreqNewContext(&entropy, &Lambda,
                                                    &iter_count,
                                                    NRrecord->mat_b,
                                                    row_probs, col_probs);
            if (status < 0) {
                return status;
            } else if (status > 0) {
                /* Failed to compute the entropy; leave the entropy
                 * unconstrained */
                status = 0;
                constrain_rel_entropy = 0;
                RE_final = 0.0;
            } else {
                RE_final = entropy;
            }
        }
        break;
    case eRelEntropyOldMatrixOldContext:
        RE_final = Blast_TargetFreqEntropy(NRrecord->mat_b);
        break;
    case eUserSpecifiedRelEntropy:
        RE_final = specifiedRE;
        break;
    default:  /* I assert that we can't get here */
        fprintf(stderr, "Unknown flag for setting relative entropy"
                "in composition matrix adjustment");
        exit(1);
    }
    Blast_ApplyPseudocounts(row_probs, length1,
                            NRrecord->first_standard_freq, pseudocounts);
    Blast_ApplyPseudocounts(col_probs, length2,
                            NRrecord->second_standard_freq, pseudocounts);

    status =
        Blast_OptimizeTargetFrequencies(&NRrecord->mat_final[0][0],
                                        COMPO_NUM_TRUE_AA,
                                        &new_iterations,
                                        &NRrecord->mat_b[0][0],
                                        row_probs, col_probs,
                                        constrain_rel_entropy,
                                        RE_final,
                                        kCompoAdjustErrTolerance,
                                        kCompoAdjustIterationLimit);
    total_iterations += new_iterations;
    if (new_iterations > max_iterations)
        max_iterations = new_iterations;

    if (status == 0) {
        status = s_ScoresStdAlphabet(matrix, NRrecord->mat_final,
                                     matrixInfo->startMatrix,
                                     row_probs, col_probs,
                                     matrixInfo->ungappedLambda);
    } else if (status == -1) {
        /* out of memory */
        status = -1;
    } else {
        /* Iteration did not converge */
        fprintf(stderr, "bad probabilities from sequence 1, length %d\n",
                length1);
        for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++)
            fprintf(stderr, "%15.12f\n", row_probs[i]);
        fprintf(stderr, "bad probabilities from sequence 2, length %d\n",
                length2);
        for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++)
            fprintf(stderr, "%15.12f\n", col_probs[i]);
        fflush(stderr);
        status = 1;
    }
    return status;
}


/* Documented in composition_adjustment.h. */
int
Blast_AdjustScores(Int4 ** matrix,
                   const Blast_AminoAcidComposition * query_composition,
                   int queryLength,
                   const Blast_AminoAcidComposition * subject_composition,
                   int subjectLength,
                   const Blast_MatrixInfo * matrixInfo,
                   ECompoAdjustModes composition_adjust_mode,
                   int RE_pseudocounts,
                   Blast_CompositionWorkspace *NRrecord,
                   EMatrixAdjustRule *matrix_adjust_rule,
                   double calc_lambda(double *,int,int,double))
{
    double LambdaRatio;      /* the ratio of the corrected
                                lambda to the original lambda */
    if (query_composition->numTrueAminoAcids == 0 ||
        subject_composition->numTrueAminoAcids == 0) {
        /* Either the query or subject contains only amibiguity
           characters, most likely because the entire subject has been
           SEGed.  Compositional adjustment is meaningless. */
        return 1;
    }
    if (matrixInfo->positionBased ||
        composition_adjust_mode == eCompositionBasedStats) {
        /* Use old-style composition-based statistics unconditionally. */
        *matrix_adjust_rule =  eCompoScaleOldMatrix;
    } else {
        /* else call Yi-Kuo's code to choose mode for matrix adjustment. */

        /* The next two arrays are letter probabilities of query and
         * match in 20 letter ARND... alphabet. */
        double permutedQueryProbs[COMPO_NUM_TRUE_AA];
        double permutedMatchProbs[COMPO_NUM_TRUE_AA];

        s_GatherLetterProbs(permutedQueryProbs, query_composition->prob);
        s_GatherLetterProbs(permutedMatchProbs, subject_composition->prob);

        *matrix_adjust_rule =
            Blast_ChooseMatrixAdjustRule(queryLength, subjectLength,
                                         permutedQueryProbs,
                                         permutedMatchProbs,
                                         matrixInfo->matrixName,
                                         composition_adjust_mode);
    }  /* end else call Yi-Kuo's code to choose mode for matrix adjustment. */

    if (eCompoScaleOldMatrix == *matrix_adjust_rule) {
        return Blast_CompositionBasedStats(matrix, &LambdaRatio, matrixInfo,
                                           query_composition->prob,
                                           subject_composition->prob,
                                           calc_lambda);
    } else {
        return
            Blast_CompositionMatrixAdj(matrix,
                                       *matrix_adjust_rule,
                                       query_composition->
                                       numTrueAminoAcids,
                                       subject_composition->
                                       numTrueAminoAcids,
                                       query_composition->prob,
                                       subject_composition->prob,
                                       RE_pseudocounts,
                                       kFixedReBlosum62,
                                       NRrecord,
                                       matrixInfo);
    }
}
