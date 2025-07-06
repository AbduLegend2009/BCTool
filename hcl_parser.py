import re
import numpy as np
import pandas as pd
from typing import Dict, List, Any, Tuple


class HCLParser:
    """Enhanced HCL parser for gene expression data"""

    def __init__(self):
        self.data = {}

    def parse_file(self, file_path: str) -> Dict[str, Any]:
        """Parse HCL file and return structured data"""
        with open(file_path, 'r') as f:
            content = f.read()

        # Remove comments
        content = re.sub(r'//.*', '', content)

        # Parse the content
        self.data = self._parse_content(content)
        return self.data

    def _parse_content(self, content: str) -> Dict[str, Any]:
        """Parse HCL content into structured data"""
        result = {}

        # Parse top-level key-value pairs
        kv_pattern = r'(\w+)\s*=\s*"([^"]*)"'
        for match in re.finditer(kv_pattern, content):
            key, value = match.groups()
            result[key] = value

        # Parse blocks
        block_pattern = r'(\w+)\s+"?([^"]*)"?\s*\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}'
        for match in re.finditer(block_pattern, content):
            block_type, block_name, block_content = match.groups()

            if block_type not in result:
                result[block_type] = {}

            parsed_block = self._parse_block(block_content)
            if block_name:
                result[block_type][block_name] = parsed_block
            else:
                result[block_type] = parsed_block

        return result

    def _parse_block(self, block_content: str) -> Dict[str, Any]:
        """Parse content within a block"""
        result = {}

        # Parse key-value pairs within block
        kv_pattern = r'(\w+)\s*=\s*(.+?)(?=\n\s*\w+\s*=|\n\s*\}|$)'
        for match in re.finditer(kv_pattern, block_content, re.DOTALL):
            key, value = match.groups()
            value = value.strip().rstrip(',')

            # Try to parse as array
            if value.startswith('[') and value.endswith(']'):
                try:
                    # Parse array values
                    array_content = value[1:-1]
                    values = [float(x.strip())
                              for x in array_content.split(',') if x.strip()]
                    result[key] = values
                except ValueError:
                    # If not all numbers, keep as strings
                    values = [x.strip().strip('"')
                              for x in array_content.split(',') if x.strip()]
                    result[key] = values
            else:
                # Try to parse as number
                value = value.strip('"')
                try:
                    if '.' in value:
                        result[key] = float(value)
                    else:
                        result[key] = int(value)
                except ValueError:
                    result[key] = value

        return result

    def to_expression_matrix(self) -> pd.DataFrame:
        """Convert parsed HCL data to gene expression matrix"""
        if 'gene_data' not in self.data:
            raise ValueError("No gene_data found in HCL file")

        gene_data = self.data['gene_data']

        # Extract expression values
        expression_matrix = []
        gene_names = []

        for gene_name, gene_info in gene_data.items():
            if 'expression_values' in gene_info:
                expression_matrix.append(gene_info['expression_values'])
                gene_names.append(gene_name)

        if not expression_matrix:
            raise ValueError("No expression values found in gene data")

        # Get condition names if available
        condition_names = None
        if 'experimental_conditions' in self.data:
            conditions = self.data['experimental_conditions']
            condition_names = [conditions.get(f'condition_{i+1}', f'Condition_{i+1}')
                               for i in range(len(expression_matrix[0]))]
        else:
            condition_names = [
                f'Condition_{i+1}' for i in range(len(expression_matrix[0]))]

        # Create DataFrame
        df = pd.DataFrame(expression_matrix, index=gene_names,
                          columns=condition_names)
        return df

    def get_metadata(self) -> Dict[str, Any]:
        """Extract metadata from parsed HCL data"""
        metadata = {}

        # Copy top-level metadata
        for key in ['version', 'dataset_name']:
            if key in self.data:
                metadata[key] = self.data[key]

        # Copy metadata block if exists
        if 'metadata' in self.data:
            metadata.update(self.data['metadata'])

        # Copy matrix config if exists
        if 'matrix_config' in self.data:
            metadata['matrix_config'] = self.data['matrix_config']

        return metadata

    def get_gene_annotations(self) -> pd.DataFrame:
        """Extract gene annotation data"""
        if 'gene_data' not in self.data:
            return pd.DataFrame()

        annotations = []
        for gene_name, gene_info in self.data['gene_data'].items():
            annotation = {'gene_name': gene_name}
            for key, value in gene_info.items():
                if key != 'expression_values':
                    annotation[key] = value
            annotations.append(annotation)

        return pd.DataFrame(annotations)


def test_hcl_parser():
    """Test the HCL parser with sample data"""
    parser = HCLParser()

    try:
        # Parse the sample file
        data = parser.parse_file('sample_data.hcl')
        print("Parsed HCL data structure:")
        print(data)

        # Convert to expression matrix
        expr_matrix = parser.to_expression_matrix()
        print(f"\nExpression Matrix Shape: {expr_matrix.shape}")
        print(expr_matrix)

        # Get metadata
        metadata = parser.get_metadata()
        print(f"\nMetadata:")
        print(metadata)

        # Get gene annotations
        annotations = parser.get_gene_annotations()
        print(f"\nGene Annotations:")
        print(annotations)

    except Exception as e:
        print(f"Error testing HCL parser: {e}")


if __name__ == "__main__":
    test_hcl_parser()
