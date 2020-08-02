import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait

    names = set(people)

    for have_trait in powerset(names):
        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    print("Input::",one_gene,two_genes,have_trait)
    probability=1
    for person in people:
        name=person

        gene = 0
        trait = False
        """Set Gene Value"""
        if name in one_gene:
            gene = 1
        elif name in two_genes:
            gene = 2
        """Set Trait Value"""
        if name in have_trait:
            trait = True

        """Parent conditions"""
        if people[person]["mother"] is None or people[person]["father"] is None:
            parent_probability=calculate_probability(gene,trait)
            print(person,gene,trait,PROBS["gene"][gene],PROBS["trait"][gene][trait],parent_probability)
            probability*=parent_probability

        else:
            """Child conditional probability"""
            father = people[person]["father"]
            mother = people[person]["mother"]
            father_inherit_probability=1
            mother_inherit_probability=1
            child_probability=1
            if father in one_gene :
                father_inherit_probability=0.5*(1-PROBS["mutation"])
            elif father in two_genes:
                father_inherit_probability=1*(1-PROBS["mutation"])
            else:
                father_inherit_probability=PROBS["mutation"]

            if mother in one_gene :
                mother_inherit_probability=0.5*(1-PROBS["mutation"])
            elif mother in two_genes:
                mother_inherit_probability=1*(1-PROBS["mutation"])
            else:
                mother_inherit_probability=PROBS["mutation"]

            if gene==2:
                child_probability=mother_inherit_probability*father_inherit_probability*PROBS["trait"][gene][trait]
            elif gene==1:
                child_probability=(mother_inherit_probability*(1-father_inherit_probability)+ father_inherit_probability * (1-mother_inherit_probability))*PROBS["trait"][gene][trait]
            else:
                child_probability=(1-mother_inherit_probability)*(1-father_inherit_probability)*PROBS["trait"][gene][trait]

            print("child_probability::",person,gene,trait,father_inherit_probability,mother_inherit_probability,PROBS["trait"][gene][trait],child_probability)
            probability*=child_probability

    print("joint_probability::",probability)
    return probability


def calculate_probability(gene,trait):
    return PROBS["gene"][gene] * PROBS["trait"][gene][trait]

def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    """print("Before Update",probabilities)"""
    for person in probabilities:
        if person in one_gene:
            probabilities[person]["gene"][1]+=p
        elif person in two_genes:
            probabilities[person]["gene"][2]+=p
        else:
            probabilities[person]["gene"][0]+=p

        if person in have_trait:
            probabilities[person]["trait"][True]+=p
        else:
            probabilities[person]["trait"][False]+=p
    """print("After Update", probabilities)"""




def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    print("Pre Normalized Probabilities",probabilities)
    for person in probabilities:
        sum=0
        for gene in probabilities[person]["gene"]:
            sum+=probabilities[person]["gene"][gene]
        normalization_factor= 1/sum
        for gene in probabilities[person]["gene"]:
            final_value= probabilities[person]["gene"][gene]* normalization_factor
            probabilities[person]["gene"][gene]=final_value
        sum=0
        for trait in probabilities[person]["trait"]:
            sum += probabilities[person]["trait"][trait]
        normalization_factor = 1 / sum
        for trait in probabilities[person]["trait"]:
            final_value= probabilities[person]["trait"][trait]* normalization_factor
            probabilities[person]["trait"][trait]=final_value

    print("Normalized Probability::",probabilities)

if __name__ == "__main__":
    main()
