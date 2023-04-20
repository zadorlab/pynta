import yaml

def check_reaction_file(file_path):
# load reaction.yaml file.
    with open(file_path,'r') as file:
        reaction_yaml = yaml.load(file, Loader=yaml.FullLoader)

# task 1: count the number of reactions in the file:
    num_reactions = len(reaction_yaml)
    print("Number of reactions:", num_reactions)

# taks 2: for each reaction, check if the molecule labesl with * (e.g. *1, *2, *3 ...) are the same in both "reactant" and "product"
# task 3: for each reaction, read reactant/product. See there are duplicates in molecule labels with *. For example, *1 cannot be used twice in either reaction/product

    for i in range(num_reactions):
        reaction = reaction_yaml[i]
        reactant = reaction['reactant']
        product = reaction['product']
    
    #check for duplicates in reactant
        reactant_molecules = reactant.split('*')[1:]
        reactant_molecules = [int(m.split()[0]) for m in reactant_molecules]
        if len(reactant_molecules) != len(set(reactant_molecules)):
            print(f'Error in reaction "index {i}": duplicate molecule labels in reactant')

    #check for duplicates in product
        product_molecules = product.split('*')[1:]
        product_molecules = [int(m.split()[0]) for m in product_molecules]
        if len(product_molecules) != len(set(product_molecules)):
            print(f'Error in reaction "index {i}": duplicate molecule labels in product')

    #check if the same molecule labels are used in reactant and product for each reaction 
        if set(reactant_molecules) != set(product_molecules):
            print(f'Warning in reaction "index {i}": molecule labels in reactant and product do not match')
    #task 4: if there are no errors, good to go!
        else:
            print(f'==== Congrats! the format for your reaction in "index {i}" looks good ====')

    return   
