#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

void get_individual_mendelian_probability_2n(int* offspring_genotype, int* sire_genotypes, int* dam_genotypes,
  int* number_sire, int* number_dam, int* number_marker, int* number_variant,
  double* allele_frequency,
  double* all_probabilities, int* all_mismatches)
{
  // individual_genotype = génotype de l'individu à assigner (vecteur de 2 * nombre de marqueurs)
  // dam_genotypes = génotypes de tous les ères possible (vecteur de 2 * nombre de marqueurs * nombre de ères)
  // dam_genotypes = génotypes de toutes les ères possible (vecteur de 2 * nombre de marqueurs * nombre de ères)

  // Pour tous les marqueurs de individual_genotype (boucle for de 0 à 2 * nombre de marqueurs de 2 en 2)
  //    génotype de l'individu au marqueur
  //    créer la table de probabilité en fonction du génotype
  //    Pour tous les pères (boucle for de 0 à nombre de ères)
  //        génotype du père au marqueur
  //        Pour toutes les ères (boucle for de 0 à nombre de ères)
  //            génotype de la mère au marqueur
  //            position dans la table
  //            probabilité mendeliennne pour le marqueur

  // VARIABLES
  int s, d, m;
  int probability_position = 0;

  // fonction
  for (s = 0; s < *number_sire; s++) { // Pour tous les pères
    // On récupère le génotype du père

    int j, pos_sire = 0;
    int initial_position = 2 * *number_marker * s;
    int final_position = initial_position + (2 * *number_marker);
    int current_sire[2 * *number_marker];

    for (j = initial_position; j < final_position; j++) {
      current_sire[pos_sire] = sire_genotypes[j];
      pos_sire++;
    }

    for (d = 0; d < *number_dam; d++) { // Pour toutes les mères
      // On récuère le génotype de la mère
      int k, pos_dam = 0;
      int initial_position = 2 * *number_marker * d;
      int final_position = initial_position + (2 * *number_marker);
      int current_dam[2 * *number_marker];

      for (k = initial_position; k < final_position; k++) {
        current_dam[pos_dam] = dam_genotypes[k];
        pos_dam++;
      }

      // Vraisemblance de chaque marqueur
      double markers_likelihood[*number_marker];
      //Mismatch de chaque marqueur
      int markers_mismatch[*number_marker];
      int marker_pos = 0;

      for (m = 0; m < 2 * *number_marker; m+=2) { // Pour chaque marqueur
        int current_individual_marker[2]  = {offspring_genotype[m], offspring_genotype[m+1]};
        int current_sire_marker[2]        = {current_sire[m], current_sire[m+1]};
        int current_dam_marker[2]         = {current_dam[m], current_dam[m+1]};

        double e = 0.01;
        int reference_allele, alternative_allele;
        double fA = 0.01, fB = 0.01;

        // Créer la table de probabilité si l'individu est homozygote
        if ((current_individual_marker[0] == current_individual_marker[1]) & (current_individual_marker[0] != 0)) {
          reference_allele =  current_individual_marker[0];
          fA = allele_frequency[(*number_variant * marker_pos) + reference_allele - 1];
        } else if (current_individual_marker[0] != current_individual_marker[1]) {
          reference_allele    = current_individual_marker[0];
          alternative_allele  = current_individual_marker[1];
          fA = allele_frequency[(*number_variant * marker_pos) + reference_allele - 1];
          fB = allele_frequency[(*number_variant * marker_pos) + alternative_allele - 1];
        }

        double homozygous_probability_table[4][4] = {{1, 0.5, e, fA},
                                                     {0.5, 0.25, e, 0.5*fA},
                                                     {e, e, e, e},
                                                     {fA, 0.5*fA, e, fA*fA}};
        double heterozygous_probability_table[7][7] = {{e,   0.5,   1,   e,    0.5,    e, fA},
                                                       {0.5, 0.5,   0.5, 0.25, 0.25,   e, 0.5*(fA+fB)},
                                                       {1,   0.5,   e,   0.5,  e,      e, fA},
                                                       {e,   0.25,  0.5, e,    0.25,   e, 0.5*fB},
                                                       {0.5, 0.25,  e,   0.25, e,      e, 0.5*fA},
                                                       {e,   e,     e,   e,    e,      e, e},
                                                       {fB,  fA+fB, fA,  fB,   0.5*fB, e, 2*fA*fB}};

        int homozygous_mismatch_table[4][4] = {{0, 0, 1, 0},
                                               {0, 0, 1, 0},
                                               {1, 1, 2, 1},
                                               {0, 0, 1, 0}};
        int heterozygous_mismatch_table[7][7] = {{1, 0, 0, 1, 0, 1, 0},
                                                 {0, 0, 0, 0, 0, 1, 0},
                                                 {0, 0, 1, 0, 1, 1, 0},
                                                 {1, 0, 0, 1, 0, 1, 0},
                                                 {0, 0, 1, 0, 1, 1, 0},
                                                 {1, 1, 1, 1, 1, 2, 1},
                                                 {0, 0, 0, 0, 0, 1, 0}};

        // Trouver les positions dans les tables en fonction des génotypes
        int sire_position, dam_position;

        // Position du père
        if ((current_individual_marker[0] == current_individual_marker[1]) & (current_individual_marker[0] == 0)) {
          // Si l'individu est NA/NA (0 0)
          sire_position = 0;
        } else if ((current_individual_marker[0] == current_individual_marker[1]) & (current_individual_marker[0] != 0)) {
          // Si l'individu est A/A
          if ((current_sire_marker[0] == current_sire_marker[1]) & (current_sire_marker[0] == 0)) {
            // père NA/NA
            sire_position = 4;
          } else if ((current_sire_marker[0] == current_sire_marker[1]) & (current_sire_marker[0] != 0)) {
            // père homozygote
            if (current_sire_marker[0] == current_individual_marker[0]) {
              // père A/A
              sire_position = 1;
            } else {
              // père C/C
              sire_position = 3;
            }
          } else {
            // père hétérozygote
            if (((current_sire_marker[0] != current_individual_marker[0]) & (current_sire_marker[1] != current_individual_marker[0])) & ((current_sire_marker[0] != current_individual_marker[1]) & (current_sire_marker[1] != current_individual_marker[1]))) {
              // père C/C
              sire_position = 3;
            } else {
              // père A/C
              sire_position = 2;
            }
          }
        } else {
          // individu hétérozygote
          if (current_sire_marker[0] == current_sire_marker[1]) {
            // père homzygous
            if (current_sire_marker[0] == 0) {
              // père NA/NA
              sire_position = 7;
            } else if (current_sire_marker[0] == current_individual_marker[0]) {
              // père A/A
              sire_position = 1;
            } else if (current_sire_marker[0] == current_individual_marker[1]) {
              // père B/B
              sire_position = 3;
            } else {
              // père C/C
              sire_position = 6;
            }
          } else {
            // père heterozygous
            if (((current_sire_marker[0] == current_individual_marker[0]) | (current_sire_marker[1] == current_individual_marker[0])) & ((current_sire_marker[0] == current_individual_marker[1]) | (current_sire_marker[1] == current_individual_marker[1]))) {
              // père A/B
              sire_position = 2;
            } else if (((current_sire_marker[0] == current_individual_marker[0]) | (current_sire_marker[1] == current_individual_marker[0])) & ((current_sire_marker[0] != current_individual_marker[1]) | (current_sire_marker[1] != current_individual_marker[1]))) {
              // père A/C
              sire_position = 4;
            } else if (((current_sire_marker[0] != current_individual_marker[0]) | (current_sire_marker[1] != current_individual_marker[0])) & ((current_sire_marker[0] == current_individual_marker[1]) | (current_sire_marker[1] == current_individual_marker[1]))) {
              // père B/C
              sire_position = 5;
            } else {
              // ère C/C
              sire_position = 6;
            }
          }
        }

        // Position de la mère
        if ((current_individual_marker[0] == current_individual_marker[1]) & (current_individual_marker[0] == 0)) {
          // Si l'individu est NA/NA (0 0)
          dam_position  = 0;
        } else if ((current_individual_marker[0] == current_individual_marker[1]) & (current_individual_marker[0] != 0)) {
          // Si l'individu est A/A
          if ((current_dam_marker[0] == current_dam_marker[1]) & (current_dam_marker[0] == 0)) {
            // mère NA/NA
            dam_position = 4;
          } else if ((current_dam_marker[0] == current_dam_marker[1]) & (current_dam_marker[0] != 0)) {
            // mère homozygote
            if (current_dam_marker[0] == current_individual_marker[0]) {
              // mère A/A
              dam_position = 1;
            } else {
              // mère C/C
              dam_position = 3;
            }
          } else {
            // mère hétérozygote
            if (((current_dam_marker[0] != current_individual_marker[0]) & (current_dam_marker[1] != current_individual_marker[0])) & ((current_dam_marker[0] != current_individual_marker[1]) & (current_dam_marker[1] != current_individual_marker[1]))) {
              // mère C/C
              dam_position = 3;
            } else {
              // mère A/C
              dam_position = 2;
            }
          }
        } else {
          // individu hétérozygote
          if (current_dam_marker[0] == current_dam_marker[1]) {
            // mère homzygous
            if (current_dam_marker[0] == 0) {
              // mère NA/NA
              dam_position = 7;
            } else if (current_dam_marker[0] == current_individual_marker[0]) {
              // mère A/A
              dam_position = 1;
            } else if (current_dam_marker[0] == current_individual_marker[1]) {
              // mère B/B
              dam_position = 3;
            } else {
              // mère C/C
              dam_position = 6;
            }
          } else {
            // mère heterozygous
            if (((current_dam_marker[0] == current_individual_marker[0]) | (current_dam_marker[1] == current_individual_marker[0])) & ((current_dam_marker[0] == current_individual_marker[1]) | (current_dam_marker[1] == current_individual_marker[1]))) {
              // mère A/B
              dam_position = 2;
            } else if (((current_dam_marker[0] == current_individual_marker[0]) | (current_dam_marker[1] == current_individual_marker[0])) & ((current_dam_marker[0] != current_individual_marker[1]) | (current_dam_marker[1] != current_individual_marker[1]))) {
              // mère A/C
              dam_position = 4;
            } else if (((current_dam_marker[0] != current_individual_marker[0]) | (current_dam_marker[1] != current_individual_marker[0])) & ((current_dam_marker[0] == current_individual_marker[1]) | (current_dam_marker[1] == current_individual_marker[1]))) {
              // mère B/C
              dam_position = 5;
            } else {
              // mère C/C
              dam_position = 6;
            }
          }
        }

        // Récupération de la Vraisemblance dans les tables
        if ((current_individual_marker[0] == current_individual_marker[1]) & (current_individual_marker[0] == 0)) {
          // individu NA/NA
          markers_likelihood[marker_pos] = 1.0;
        } else if ((current_individual_marker[0] == current_individual_marker[1]) & (current_individual_marker[0] != 0)) {
          // individu homozygote
          markers_likelihood[marker_pos] = homozygous_probability_table[sire_position - 1][dam_position - 1];
        } else {
          // individu hétérozygote
          markers_likelihood[marker_pos] = heterozygous_probability_table[sire_position - 1][dam_position - 1];
        }

        // Récupération des mismatches dans les tables
        if ((current_individual_marker[0] == current_individual_marker[1]) & (current_individual_marker[0] == 0)) {
          // individu NA/NA
          markers_mismatch[marker_pos] = 0;
        } else if ((current_individual_marker[0] == current_individual_marker[1]) & (current_individual_marker[0] != 0)) {
          // individu homozygote
          markers_mismatch[marker_pos] = homozygous_mismatch_table[sire_position - 1][dam_position - 1];
        } else {
          // individu hétérozygote
          markers_mismatch[marker_pos] = heterozygous_mismatch_table[sire_position - 1][dam_position - 1];
        }

        /*
        printf("ind : %d - %d | sire : %d - %d | dam : %d - %d | P = %f",
        current_individual_marker[0], current_individual_marker[1],
        current_sire_marker[0], current_sire_marker[1],
        current_dam_marker[0], current_dam_marker[1],
        markers_likelihood[marker_pos]);
        printf("\n");
        */


        marker_pos += 1;
      }

      double V = 0.0;
      int R = 0;
      for (m = 0; m < *number_marker; m++) {
        V = V + log(markers_likelihood[m]);
        R = R + markers_mismatch[m];
      }

      double P = exp(V / *number_marker);
      all_probabilities[probability_position] = P;
      all_mismatches[probability_position] = R;

      probability_position++;
    }
  }
}

void get_individual_mendelian_probability_3n(int* offspring_genotype, int* sire_genotypes, int* dam_genotypes,
  int* number_sire, int* number_dam, int* number_marker, int* number_variant,
  double* allele_frequency,
  double* all_probabilities, int* all_mismatches,
  double* t_recom)
{
  // individual_genotype = génotype de l'individu à assigner (vecteur de 2 * nombre de marqueurs)
  // dam_genotypes = génotypes de tous les ères possible (vecteur de 2 * nombre de marqueurs * nombre de ères)
  // dam_genotypes = génotypes de toutes les ères possible (vecteur de 2 * nombre de marqueurs * nombre de ères)

  // Pour tous les marqueurs de individual_genotype (boucle for de 0 à 2 * nombre de marqueurs de 2 en 2)
  //    génotype de l'individu au marqueur
  //    créer la table de probabilité en fonction du génotype
  //    Pour tous les pères (boucle for de 0 à nombre de ères)
  //        génotype du père au marqueur
  //        Pour toutes les ères (boucle for de 0 à nombre de ères)
  //            génotype de la mère au marqueur
  //            position dans la table
  //            probabilité mendeliennne pour le marqueur

  // VARIABLES
  int s, d, m;
  int probability_position = 0;
  double r = *t_recom;

  for (s = 0; s < *number_sire; s++) { // Pour tous les pères
    // On récupère le génotype du père

    int j, pos_sire = 0;
    int initial_position = 2 * *number_marker * s;
    int final_position = initial_position + (2 * *number_marker);
    int current_sire[2 * *number_marker];

    for (j = initial_position; j < final_position; j++) {
      current_sire[pos_sire] = sire_genotypes[j];
      pos_sire++;
    }

    for (d = 0; d < *number_dam; d++) { // Pour toutes les mères
      // On récuère le génotype de la mère
      int k, pos_dam = 0;
      int initial_position = 2 * *number_marker * d;
      int final_position = initial_position + (2 * *number_marker);
      int current_dam[2 * *number_marker];

      for (k = initial_position; k < final_position; k++) {
        current_dam[pos_dam] = dam_genotypes[k];
        pos_dam++;
      }

      // Vraisemblance de chaque marqueur
      double markers_likelihood[*number_marker];
      //Mismatch de chaque marqueur
      int markers_mismatch[*number_marker];
      int marker_pos = 0;

      // Parental position
      int diploid_pos = 0;

      int allele_A = 0, allele_B = 0, allele_C = 0;

      for (m = 0; m < 3 * *number_marker; m+=3) { // Pour chaque marqueur
        int current_individual_marker[3]  = {offspring_genotype[m], offspring_genotype[m+1], offspring_genotype[m+2]};
        int current_sire_marker[2]        = {current_sire[diploid_pos], current_sire[diploid_pos+1]};
        int current_dam_marker[2]         = {current_dam[diploid_pos], current_dam[diploid_pos+1]};
        diploid_pos += 2;

        double e = 0.01; // genotyping error rate
        // double r = 0.5; // recombination rate
        int reference_allele, alternative_allele, off_allele;
        double fA = 0.01, fB = 0.01, fC = 0.01;

        // Créer la table de probabilité si l'individu est homozygote (A/A/A)
        if ((current_individual_marker[0] == current_individual_marker[1]) & (current_individual_marker[1] == current_individual_marker[2]) & (current_individual_marker[0] != 0)) {
          reference_allele =  current_individual_marker[0];
          fA = allele_frequency[(*number_variant * marker_pos) + reference_allele - 1];
        } else if ((current_individual_marker[0] != current_individual_marker[1]) & (current_individual_marker[1] == current_individual_marker[2])) {
          // individual is A/B/B
          reference_allele    = current_individual_marker[0];
          alternative_allele  = current_individual_marker[1];
          fA = allele_frequency[(*number_variant * marker_pos) + reference_allele - 1];
          fB = allele_frequency[(*number_variant * marker_pos) + alternative_allele - 1];
        } else if ((current_individual_marker[0] == current_individual_marker[1]) & (current_individual_marker[1] != current_individual_marker[2])) {
          // individual is A/A/B
          reference_allele    = current_individual_marker[2];
          alternative_allele  = current_individual_marker[0];
          fA = allele_frequency[(*number_variant * marker_pos) + reference_allele - 1];
          fB = allele_frequency[(*number_variant * marker_pos) + alternative_allele - 1];
        } else {
          // individual is A/B/C
          reference_allele    = current_individual_marker[0];
          alternative_allele  = current_individual_marker[1];
          off_allele          = current_individual_marker[2];
          fA = allele_frequency[(*number_variant * marker_pos) + reference_allele - 1];
          fB = allele_frequency[(*number_variant * marker_pos) + alternative_allele - 1];
          fC = allele_frequency[(*number_variant * marker_pos) + off_allele - 1];
        }

        // AA - AD - DD - missing
        double AAA_probability_table[4][4] = {
          {1,   0.5*(1-r),    e, r*fA*fA+(1-r)*fA},
          {0.5, 0.25*(1-r),   e, 0.5*(r*fA*fA + (1-r)*fA)},
          {e,   e,        e, e},
          {fA,  0.5*(1-r)*fA,  e, r*fA*fA*fA+(1-r)*fA*fA}
        };

        // AA - AB - BB - AD - BD - DD - missing
        double ABB_probability_table[7][7] = {
          {e, 0.5*(1-r),            1,   e, 0.5*(1-r),      e, (1-r)*fA*fB+fB*fB+(1-r)*fB*fC},
          {e, 0.25*(1-r) + 0.5*r, 0.5, e, 0.25*(1-r),  e, (1-0.5*(1-r))*fA*fB+0.5*fB*fB+0.5*(1-r)*fB*fC},
          {e, r,                    e,   e, e,              e, 2*fA*fB*r},
          {e, 0.25*(1-r),           0.5, e, 0.25*(1-r),     e, 0.5*(1-r)*fA*fB+0.5*(1-r)*fB*fC+0.5*fB*fB},
          {e, 0.5*r,                e,   e, e,              e, fA*fB*r},
          {e, e,                    e,   e, e,              e, e},
          {e, 0.5*(1-r)*fA+r*fB,  fA,  e, 0.5*fA*(1-r),   e, fA*((1-r)*fA*fB+fB*fB+(1-r)*fB*fC)+fB*fB*fA*2*r}
        };

        // AA	AB	BB	AC	BC	CC	AD	BD	CD	DD	missing
        double ABC_probability_table[11][11] = {
          {e, e,      e, e,       r,      e, e, e, e, e, 2*fB*fC*r},
          {e, e,      e, 0.5*r,   0.5*r,  e, e, e, e, e, fB*fC*r+fA*fC*r},
          {e, e,      e, r,       e,      e, e, e, e, e, 2*fA*fC*r},
          {e, 0.5*r,  e, e ,      0.5*r,  e, e, e, e, e, fA*fB*r+fB*fC*r},
          {e, 0.5*r,  e, 0.5*r,   e,      e, e, e, e, e, fA*fB*r+fA*fC*r},
          {e, r,      e, e,       e,      e, e, e, e, e, 2*fA*fB*r},
          {e, e,      e, e,       0.5*r,  e, e, e, e, e, fB*fC*r},
          {e, e,      e, 0.5*r,   e,      e, e, e, e, e, fA*fC*r},
          {e, 0.5*r,  e, e,       e,      e, e, e, e, e, fB*fA*r},
          {e, e,      e, e,       e,      e, e, e, e, e, e},
          {e, fC*0.5*r,   e, fB*0.5*r,    fA*0.5*r,   e, e, e, e, e, 6*fA*fB*fC*r}
        };

        // AA - AD - DD - missing
        int AAA_mismatch_table[4][4] = {
          {0, 0, 1, 0},
          {0, 0, 1, 0},
          {1, 1, 2, 1},
          {0, 0, 1, 0}
       };

        // AA - AB - BB - AD - BD - DD - missing
        int ABB_mismatch_table[7][7] = {
          {1, 0, 0, 1, 0, 1, 0},
          {1, 0, 0, 1, 0, 1, 0},
          {1, 0, 1, 1, 1, 1, 0},
          {1, 0, 0, 1, 0, 1, 0},
          {1, 0, 1, 1, 1, 1, 0},
          {2, 1, 1, 2, 1, 2, 1},
          {1, 0, 0, 1, 0, 1, 0}
        };

        // AA	AB	BB	AC	BC	CC	AD	BD	CD	DD	missing
        int ABC_mismatch_table[11][11] = {
          {1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0},
          {1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0},
          {1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0},
          {1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0},
          {1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0},
          {1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0},
          {1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0},
          {1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0},
          {1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0},
          {2, 1, 2, 1, 1, 2, 2, 2, 2, 2, 1},
          {1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0}
        };

        // Trouver les positions dans les tables en fonction des génotypes
        int sire_position, dam_position;

        // Obtenir les trois allèles de l'individu
        allele_A = 0;
        allele_B = 0;
        allele_C = 0;

        if (current_individual_marker[0] == 0) {
          // individual is "NA/NA/NA"
          allele_A = 0;
          allele_B = 0;
          allele_C = 0;
        } else if ((current_individual_marker[0] == current_individual_marker[1]) & (current_individual_marker[1] == current_individual_marker[2]) & (current_individual_marker[0] != 0)) {
          // individual is homozygous
          allele_A = current_individual_marker[0];
        } else if ((current_individual_marker[0] != current_individual_marker[1]) & (current_individual_marker[1] == current_individual_marker[2])) {
          // individual is "A/B/B"
          allele_A = current_individual_marker[0];
          allele_B = current_individual_marker[1];
        } else if ((current_individual_marker[0] == current_individual_marker[1]) & (current_individual_marker[1] != current_individual_marker[2])) {
          // individual is A/A/B
          allele_A = current_individual_marker[2];
          allele_B = current_individual_marker[0];
        } else {
          // individual is "A/B/C"
          allele_A = current_individual_marker[0];
          allele_B = current_individual_marker[1];
          allele_C = current_individual_marker[2];
        }

        // Obtenir les deux allèles du père
        bool sire_is_AA = false;
        bool sire_is_AB = false;
        bool sire_is_AC = false;
        bool sire_is_AD = false;
        bool sire_is_BB = false;
        bool sire_is_BC = false;
        bool sire_is_BD = false;
        bool sire_is_CC = false;
        bool sire_is_CD = false;
        bool sire_is_DD = false;

        bool sire_is_NA = false;

        if (current_sire_marker[0] == 0) {
          // sire is "NA/NA"
          sire_is_NA = true;
        } else if ((current_sire_marker[0] == current_sire_marker[1]) & (current_sire_marker[0] != 0)) {
          // sire is homozygous
          if ((allele_A == 0)) {
            // individual is "NA/NA/NA"
          } else if ((allele_A != 0) & (allele_B == 0) & (allele_C == 0)) {
            // individual is "A/A/A"
            if (current_sire_marker[0] == 0) {
              sire_is_NA = true;
            } else if (current_sire_marker[0] == allele_A) {
              sire_is_AA = true;
            } else {
              sire_is_DD = true;
            }
          } else if ((allele_A != 0) & (allele_B != 0) & (allele_C == 0)) {
            // individual is "A/B/B"
            if (current_sire_marker[0] == 0) {
              sire_is_NA = true;
            } else if (current_sire_marker[0] == allele_A) {
              sire_is_AA = true;
            } else if (current_sire_marker[0] == allele_B) {
              sire_is_BB = true;
            } else {
              sire_is_DD = true;
            }
          } else {
            // individual is "A/B/C"
            if (current_sire_marker[0] == 0) {
              sire_is_NA = true;
            } else if (current_sire_marker[0] == allele_A) {
              sire_is_AA = true;
            } else if (current_sire_marker[0] == allele_B) {
              sire_is_BB = true;
            } else if (current_sire_marker[0] == allele_C) {
              sire_is_CC = true;
            } else {
              sire_is_DD = true;
            }
          }
        } else {
          // sire is heterozygous
          if ((allele_A != 0) & (allele_B == 0) & (allele_C == 0)) {
            // individual is "A/A/A"
            if ((current_sire_marker[0] == allele_A) | (current_sire_marker[1] == allele_A)) {
              sire_is_AD = true;
            }
          } else if ((allele_A != 0) & (allele_B != 0) & (allele_C == 0)) {
            // individual is "A/B/B"
            if ((current_sire_marker[0] == allele_A) & (current_sire_marker[1] == allele_B)) {
              sire_is_AB = true;
            } else if ((current_sire_marker[0] == allele_B) & (current_sire_marker[1] == allele_A)) {
              sire_is_AB = true;
            } else if ((current_sire_marker[0] == allele_A) & (current_sire_marker[1] != allele_B)) {
              sire_is_AD = true;
            } else if ((current_sire_marker[0] != allele_A) & (current_sire_marker[1] == allele_B)) {
              sire_is_AD = true;
            } else if ((current_sire_marker[0] == allele_B) & (current_sire_marker[1] != allele_A)) {
              sire_is_BD = true;
            } else if ((current_sire_marker[0] != allele_B) & (current_sire_marker[1] == allele_A)) {
              sire_is_BD = true;
            } else {
              sire_is_DD = true;
            }
          } else {
            if ((current_sire_marker[0] == allele_A) & (current_sire_marker[1] == allele_B) & (allele_A != 0) & (allele_B != 0)) {
              sire_is_AB = true;
            } else if ((current_sire_marker[0] == allele_A) & (current_sire_marker[1] == allele_C) & (allele_A != 0) & (allele_C != 0)) {
              sire_is_AC = true;
            } else if ((current_sire_marker[0] == allele_A) & (current_sire_marker[1] != allele_B) & (current_sire_marker[1] != allele_C) & (allele_A != 0) & (allele_B != 0) & (allele_C != 0)) {
              sire_is_AD = true;
            } else if ((current_sire[0] == allele_B) & (current_sire_marker[1] == allele_A) & (allele_B != 0) & (allele_A != 0)) {
              sire_is_AB = true;
            } else if ((current_sire[0] == allele_B) & (current_sire_marker[1] == allele_C) & (allele_B != 0) & (allele_C != 0)) {
              sire_is_BC = true;
            } else if ((current_sire[0] == allele_B) & (current_sire_marker[1] != allele_A) & (current_sire_marker[1] != allele_C) & (allele_B != 0) & (allele_A != 0) & (allele_C != 0)) {
              sire_is_BD = true;
            } else if ((current_sire[0] == allele_C) & (current_sire[1] == allele_A) & (allele_C != 0) & (allele_A != 0)) {
              sire_is_AC = true;
            } else if ((current_sire[0] == allele_C) & (current_sire[1] == allele_B) & (allele_C != 0) & (allele_B != 0)) {
              sire_is_BC = true;
            } else if ((current_sire[0] == allele_C) & (current_sire[1] != allele_A) & (current_sire[1] != allele_B) & (allele_C != 0) & (allele_A != 0) & (allele_B != 0)) {
              sire_is_CD = true;
            } else if ((current_sire[0] != allele_A) & (current_sire[0] != allele_B) & (current_sire[0] != allele_C) & (current_sire[1] != allele_A) & (current_sire[1] != allele_B) & (current_sire[1] != allele_C) & (allele_C != 0) & (allele_A != 0) & (allele_B != 0)) {
              sire_is_DD = true;
            }
          }
        }

        // Obtenir les deux allèles du mère
        bool dam_is_AA = false;
        bool dam_is_AB = false;
        bool dam_is_AC = false;
        bool dam_is_AD = false;
        bool dam_is_BB = false;
        bool dam_is_BC = false;
        bool dam_is_BD = false;
        bool dam_is_CC = false;
        bool dam_is_CD = false;
        bool dam_is_DD = false;

        bool dam_is_NA = false;

        if (current_dam_marker[0] == 0) {
          // dam is "NA/NA"
          dam_is_NA = true;
        } else if ((current_dam_marker[0] == current_dam_marker[1]) & (current_dam_marker[0] != 0)) {
          // dam is homozygous
          if ((allele_A == 0)) {
            // individual is "NA/NA/NA"
          } else if ((allele_A != 0) & (allele_B == 0) & (allele_C == 0)) {
            // individual is "A/A/A"
            if (current_dam_marker[0] == 0) {
              dam_is_NA = true;
            } else if (current_dam_marker[0] == allele_A) {
              dam_is_AA = true;
            } else {
              dam_is_DD = true;
            }
          } else if ((allele_A != 0) & (allele_B != 0) & (allele_C == 0)) {
            // individual is "A/B/B"
            if (current_dam_marker[0] == 0) {
              dam_is_NA = true;
            } else if (current_dam_marker[0] == allele_A) {
              dam_is_AA = true;
            } else if (current_dam_marker[0] == allele_B) {
              dam_is_BB = true;
            } else {
              dam_is_DD = true;
            }
          } else {
            // individual is "A/B/C"
            if (current_dam_marker[0] == 0) {
              dam_is_NA = true;
            } else if (current_dam_marker[0] == allele_A) {
              dam_is_AA = true;
            } else if (current_dam_marker[0] == allele_B) {
              dam_is_BB = true;
            } else if (current_dam_marker[0] == allele_C) {
              dam_is_CC = true;
            } else {
              dam_is_DD = true;
            }
          }
        } else {
          // dam is heterozygous
          if ((allele_A != 0) & (allele_B == 0) & (allele_C == 0)) {
            // individual is "A/A/A"
            if ((current_dam_marker[0] == allele_A) | (current_dam_marker[1] == allele_A)) {
              dam_is_AD = true;
            }
          } else if ((allele_A != 0) & (allele_B != 0) & (allele_C == 0)) {
            // individual is "A/B/B"
            if ((current_dam_marker[0] == allele_A) & (current_dam_marker[1] == allele_B)) {
              dam_is_AB = true;
            } else if ((current_dam_marker[0] == allele_B) & (current_dam_marker[1] == allele_A)) {
              dam_is_AB = true;
            } else if ((current_dam_marker[0] == allele_A) & (current_dam_marker[1] != allele_B)) {
              dam_is_AD = true;
            } else if ((current_dam_marker[0] != allele_A) & (current_dam_marker[1] == allele_B)) {
              dam_is_AD = true;
            } else if ((current_dam_marker[0] == allele_B) & (current_dam_marker[1] != allele_A)) {
              dam_is_BD = true;
            } else if ((current_dam_marker[0] != allele_B) & (current_dam_marker[1] == allele_A)) {
              dam_is_BD = true;
            } else {
              dam_is_DD = true;
            }
          } else {
            if ((current_dam_marker[0] == allele_A) & (current_dam_marker[1] == allele_B) & (allele_A != 0) & (allele_B != 0)) {
              dam_is_AB = true;
            } else if ((current_dam_marker[0] == allele_A) & (current_dam_marker[1] == allele_C) & (allele_A != 0) & (allele_C != 0)) {
              dam_is_AC = true;
            } else if ((current_dam_marker[0] == allele_A) & (current_dam_marker[1] != allele_B) & (current_dam_marker[1] != allele_C) & (allele_A != 0) & (allele_B != 0) & (allele_C != 0)) {
              dam_is_AD = true;
            } else if ((current_dam[0] == allele_B) & (current_dam_marker[1] == allele_A) & (allele_B != 0) & (allele_A != 0)) {
              dam_is_AB = true;
            } else if ((current_dam[0] == allele_B) & (current_dam_marker[1] == allele_C) & (allele_B != 0) & (allele_C != 0)) {
              dam_is_BC = true;
            } else if ((current_dam[0] == allele_B) & (current_dam_marker[1] != allele_A) & (current_dam_marker[1] != allele_C) & (allele_B != 0) & (allele_A != 0) & (allele_C != 0)) {
              dam_is_BD = true;
            } else if ((current_dam[0] == allele_C) & (current_dam[1] == allele_A) & (allele_C != 0) & (allele_A != 0)) {
              dam_is_AC = true;
            } else if ((current_dam[0] == allele_C) & (current_dam[1] == allele_B) & (allele_C != 0) & (allele_B != 0)) {
              dam_is_BC = true;
            } else if ((current_dam[0] == allele_C) & (current_dam[1] != allele_A) & (current_dam[1] != allele_B) & (allele_C != 0) & (allele_A != 0) & (allele_B != 0)) {
              dam_is_CD = true;
            } else if ((current_dam[0] != allele_A) & (current_dam[0] != allele_B) & (current_dam[0] != allele_C) & (current_dam[1] != allele_A) & (current_dam[1] != allele_B) & (current_dam[1] != allele_C) & (allele_C != 0) & (allele_A != 0) & (allele_B != 0)) {
              dam_is_DD = true;
            }
          }
        }


        // Obtenir les positions
        if ((allele_A == 0) & (allele_B == 0) & (allele_C == 0)) {
          // individual is NA/NA
          sire_position = 0;
          dam_position = 0;
        } else if ((allele_A != 0) & (allele_B == 0) & (allele_C == 0)) {
          // individual is "A/A/A"
          if (sire_is_AA) {
            sire_position = 1;
          } else if (sire_is_AD) {
            sire_position = 2;
          } else if (sire_is_DD) {
            sire_position = 3;
          } else if (sire_is_NA) {
            sire_position = 4;
          } else {
            sire_position = 4;
          }

          if (dam_is_AA) {
            dam_position = 1;
          } else if (dam_is_AD) {
            dam_position = 2;
          } else if (dam_is_DD) {
            dam_position = 3;
          } else if (dam_is_NA) {
            dam_position = 4;
          } else {
            dam_position = 4;
          }
        } else if ((allele_A != 0) & (allele_B != 0) & (allele_C == 0)) {
          // individual is "A/B/B"
          if (sire_is_AA) {
            sire_position = 1;
          } else if (sire_is_AB) {
            sire_position = 2;
          } else if (sire_is_BB) {
            sire_position = 3;
          } else if (sire_is_AD) {
            sire_position = 4;
          } else if (sire_is_BD) {
            sire_position = 5;
          } else if (sire_is_DD) {
            sire_position = 6;
          } else if (sire_is_NA) {
            sire_position = 7;
          } else {
            sire_position = 7;
          }

          if (dam_is_AA) {
            dam_position = 1;
          } else if (dam_is_AB) {
            dam_position = 2;
          } else if (dam_is_BB) {
            dam_position = 3;
          } else if (dam_is_AD) {
            dam_position = 4;
          } else if (dam_is_BD) {
            dam_position = 5;
          } else if (dam_is_DD) {
            dam_position = 6;
          } else if (dam_is_NA) {
            dam_position = 7;
          } else {
            dam_position = 7;
          }
        } else {
          // individual is "A/B/C"
          if (sire_is_AA) {
            sire_position = 1;
          } else if (sire_is_AB) {
            sire_position = 2;
          } else if (sire_is_BB) {
            sire_position = 3;
          } else if (sire_is_AC) {
            sire_position = 4;
          } else if (sire_is_BC) {
            sire_position = 5;
          } else if (sire_is_CC) {
            sire_position = 6;
          } else if (sire_is_AD) {
            sire_position = 7;
          } else if (sire_is_BD) {
            sire_position = 8;
          } else if (sire_is_CD) {
            sire_position = 9;
          } else if (sire_is_DD) {
            sire_position = 10;
          } else if (sire_is_NA) {
            sire_position = 11;
          } else {
            sire_position = 11;
          }

          if (dam_is_AA) {
            dam_position = 1;
          } else if (dam_is_AB) {
            dam_position = 2;
          } else if (dam_is_BB) {
            dam_position = 3;
          } else if (dam_is_AC) {
            dam_position = 4;
          } else if (dam_is_BC) {
            dam_position = 5;
          } else if (dam_is_CC) {
            dam_position = 6;
          } else if (dam_is_AD) {
            dam_position = 7;
          } else if (dam_is_BD) {
            dam_position = 8;
          } else if (dam_is_CD) {
            dam_position = 9;
          } else if (dam_is_DD) {
            dam_position = 10;
          } else if (dam_is_NA) {
            dam_position = 11;
          } else {
            dam_position = 11;
          }
        }

        // Obtenir les valeurs
        if ((allele_A == 0) & (allele_B == 0) & (allele_C == 0)) {
          // individual is NA/NA
          markers_likelihood[marker_pos] = 1.0;
          markers_mismatch[marker_pos] = 0;
        } else if ((allele_A != 0) & (allele_B == 0) & (allele_C == 0)) {
          // individual is "A/A/A"
          markers_likelihood[marker_pos]  = AAA_probability_table[sire_position - 1][dam_position - 1];
          markers_mismatch[marker_pos]    = AAA_mismatch_table[sire_position - 1][dam_position - 1];
        } else if ((allele_A != 0) & (allele_B != 0) & (allele_C == 0)) {
          // individual is "A/B/B"
          markers_likelihood[marker_pos]  = ABB_probability_table[sire_position - 1][dam_position - 1];
          markers_mismatch[marker_pos]    = ABB_mismatch_table[sire_position - 1][dam_position - 1];
        } else {
          // individual is "A/B/C"
          markers_likelihood[marker_pos]  = ABC_probability_table[sire_position - 1][dam_position - 1];
          markers_mismatch[marker_pos]    = ABC_mismatch_table[sire_position - 1][dam_position - 1];
        }

        /*
        printf("=================================================================\n");
        printf("marker : %d", marker_pos);
        printf("\n");
        printf("allele A : %d | allele_B : %d | allele_C : %d",
        allele_A, allele_B, allele_C);
        printf("\n");
        printf("ind : %d - %d - %d | sire : %d - %d | dam : %d - %d",
        current_individual_marker[0], current_individual_marker[1], current_individual_marker[2],
        current_sire_marker[0], current_sire_marker[1],
        current_dam_marker[0], current_dam_marker[1]);
        printf("\n");
        printf("sire AA : %d | AB : %d | BB : %d | AC : %d | BC : %d | CC : %d | AD : %d | BD : %d | BC : %d | DD : %d | NA : %d",
        sire_is_AA, sire_is_AB, sire_is_BB, sire_is_AC, sire_is_BC, sire_is_CC, sire_is_AD, sire_is_BD, sire_is_CD, sire_is_DD, sire_is_NA);
        printf("\n");
        printf("dam AA : %d | AB : %d | BB : %d | AC : %d | BC : %d | CC : %d | AD : %d | BD : %d | BC : %d | DD : %d | NA : %d",
        dam_is_AA, dam_is_AB, dam_is_BB, dam_is_AC, dam_is_BC, dam_is_CC, dam_is_AD, dam_is_BD, dam_is_CD, dam_is_DD, dam_is_NA);
        printf("\n");
        printf("sire position : %d | dam position : %d | proba : %f",
        sire_position, dam_position, markers_likelihood[marker_pos]);
        printf("\n");
        printf("fA : %f | fB : %f | fC : %f",
        fA, fB, fC);
        printf("\n");
        */


        marker_pos += 1;
      }

      double V = 0.0;
      int R = 0;
      for (m = 0; m < *number_marker; m++) {
        V = V + log(markers_likelihood[m]);
        R = R + markers_mismatch[m];
      }

      double P = exp(V / *number_marker);
      all_probabilities[probability_position] = P;
      all_mismatches[probability_position] = R;

      probability_position++;
    }
  }
}
