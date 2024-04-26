from typing import Union, Tuple
import ast
from cobra.core.gene import GPR
import cobra
import numpy as np 
import pandas as pd 

def remove_gene_from_ast(node, knockout_gene):
    if isinstance(node, ast.BoolOp):
        node.values = [remove_gene_from_ast(value, knockout_gene) for value in node.values if value is not None]
        return node if node.values else None
    elif isinstance(node, ast.Name) and node.id == knockout_gene:
        return None
    elif isinstance(node, ast.Attribute):
        if isinstance(node.value, ast.Name) and node.value.id == knockout_gene:
            return None
    return node



def _eval_gpr(
    expr: ast.Expression,
    knockouts: set,
    gene_values: dict
) -> Union[bool, Tuple[float, str]]:
    nan_genes = {key for key, value in gene_values.items() if np.isnan(value)}
    if(len(nan_genes) > 0):
        for gene in nan_genes:
            expr = remove_gene_from_ast(expr, gene)
    for gene in knockouts:
        expr = remove_gene_from_ast(expr, gene)
    
    if isinstance(expr, ast.BoolOp):
        op = expr.op
        if isinstance(op, ast.Or):
            max_value = max(_eval_gpr(i, knockouts, gene_values) for i in expr.values)
            return max_value if isinstance(max_value, tuple) else (max(max_value), "")
        elif isinstance(op, ast.And):
            min_value = min(_eval_gpr(i, knockouts, gene_values) for i in expr.values)
            if isinstance(min_value, tuple):
                return min_value[0], min_value[1]
            else:
                return min(1.0, min_value), ""
        else:
            raise TypeError(f"Unsupported operation: {op.__class__.__name__}")
    elif isinstance(expr, ast.Name):
        gene_id = expr.id.replace('@', '_')
        return gene_values.get(gene_id, 0.0), gene_id

    elif expr is None:
        return 0.0, ""
    else:
        raise TypeError(f"Unsupported operation: {repr(expr)}")


def evaluate_gpr(
    expr: str,
    knockouts: set,
    gene_values: dict
) -> Union[bool, Tuple[float, str]]:
    expression_ast = GPR.from_string(expr)
    result = _eval_gpr(expr=expression_ast.body, knockouts=knockouts, gene_values=gene_values)
    return result

expression = "56301_AT1 and 6520_AT1"
# expression = "6520_AT1"

values = {'6520_AT1': 5466.8, '56301_AT1': np.nan}
# values = {'6520_AT1': 5466.8}

# result = evaluate_gpr(expr=expression, knockouts={'56301_AT1'}, gene_values=values)
# print(result)
