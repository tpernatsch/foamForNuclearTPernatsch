# run.py
from case import build_case

# Parameters defined here
enrichments = [0.046, 0.02] 

# Run function
def run_case(
    case_name: str,
    enrichment: float
):
    """Build case with given (enrichment, case_name), export & run it."""
    case, rod_mesh = build_case(case_name=case_name, enrichment=enrichment)

    print(f"Running case: enrichment={enrichment}")
    case.clean()
    case.run()

    return case, rod_mesh


if __name__ == "__main__":
    
    for e in enrichments:
        case_name = f"cases/e{e}"

        case, rod_mesh = run_case(
            case_name=case_name,
            enrichment=e
        )