from neo4j import GraphDatabase
import csv
import argparse

uri = "bolt://localhost:7687"
username = "neo4j"
password = "password"

parser = argparse.ArgumentParser(description="Query Neo4j database and export results to CSV")
parser.add_argument("-l", required=True, help="Specify the label(s) as a start:end range.")
args = parser.parse_args()
# Parse the label range
label_range = args.l.split(':')
start_label = int(label_range[0])
end_label = int(label_range[1])

# Create a list of labels from the range
labels = [str(i) for i in range(start_label, end_label + 1)]

# Iterate through each label in BSW
for myLabel in labels:
    # Add label to selected nodes
    add_label_cypher = f"""
    MATCH (n)
    WHERE n.no_of_samples <> 1 AND n.no_of_samples <> 103 AND (n.calc = 'BSW' OR n.calc = 'standard')
    WITH n.img_name AS imgName, COLLECT(n) AS nodes
    WITH imgName, nodes, nodes[TOINTEGER(RAND() * SIZE(nodes))] AS randomNode
    SET randomNode:BSW_it{myLabel}
    """

    # Create graph projection
    graph_projection_cypher = f"""
    CALL gds.graph.project(
      'BSW_it{myLabel}',
      'BSW_it{myLabel}',
      {{
        RELATIONSHIP: {{
          orientation: 'UNDIRECTED',
          type: '*',
          properties: {{
            num_common_samples: {{
              property: 'num_common_samples'
            }}
          }}
        }}
      }}
    )
    """

    # Modularity optimization community detection
    modopt_cypher = f"""
    CALL gds.modularityOptimization.stream('BSW_it{myLabel}', {{ relationshipWeightProperty: 'num_common_samples' }})
    YIELD nodeId, communityId
    RETURN nodeId, communityId,
    gds.util.asNode(nodeId).img_name AS img_name, 
    gds.util.asNode(nodeId).tissue_region AS tissue_region,
    gds.util.asNode(nodeId).event AS event,
    gds.util.asNode(nodeId).calc as calc,
    gds.util.asNode(nodeId).leaf_status AS leaf_status,
    gds.util.asNode(nodeId).median_signal AS median_signal,
    gds.util.asNode(nodeId).molecular_data AS molecular_data,
    gds.util.asNode(nodeId).no_of_samples AS no_of_samples,
    gds.util.asNode(nodeId).quartile AS quartile,
    gds.util.asNode(nodeId).samples AS samples,
    gds.util.asNode(nodeId).subsample AS subsample
    ORDER BY communityId DESC;
    """

    csv_file_path = f"/Users/m254284/Desktop/Neo4j_Projects/Img_events/ModOpt_bycalc/Iterations/BSW/modopt_calc_it{myLabel}.csv"

    # Open the CSV file for writing
    with open(csv_file_path, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        # Write header row
        csv_writer.writerow(['node_id','community','img_name', 'tissue_region','event','calc','leaf_status',
                             'median_signal', 'molecular_data','no_of_samples', 'quartile', 'samples', 'subsample'])

        with GraphDatabase.driver(uri, auth=(username, password)) as driver:
            with driver.session() as session:
                # Add label to selected nodes
                session.run(add_label_cypher)

                # Create graph projection
                session.run(graph_projection_cypher)

                # Run modularity optimization community detection
                result = session.run(modopt_cypher)

                # Write the result to the CSV file
                for record in result:
                    node_id = record['nodeId']
                    community = record['communityId']
                    img_name = record['img_name']
                    tissue_region = record['tissue_region']
                    event = record['event']
                    calc = record['calc']
                    leaf_status = record['leaf_status']
                    median_signal = record['median_signal']
                    molecular_data = record['molecular_data']
                    no_of_samples = record['no_of_samples']
                    quartile = record['quartile']
                    samples = record['samples']
                    subsample = record['subsample']

                    csv_writer.writerow([node_id, community, img_name, tissue_region, event, calc, leaf_status,
                                     median_signal, molecular_data, no_of_samples, quartile, samples, subsample])

    # Close the Neo4j driver
    driver.close()

    # Print a message indicating completion for each label
    print(f"Finished all communities in BSW_it{myLabel}.")


# ----------- QWS
# Iterate through each label in BSW
for myLabel in labels:
    # Add label to selected nodes
    add_label_cypher = f"""
    MATCH (n)
    WHERE n.no_of_samples <> 1 AND n.no_of_samples <> 103 AND (n.calc = 'QWS' OR n.calc = 'standard')
    WITH n.img_name AS imgName, COLLECT(n) AS nodes
    WITH imgName, nodes, nodes[TOINTEGER(RAND() * SIZE(nodes))] AS randomNode
    SET randomNode:QWS_it{myLabel}
    """

    # Create graph projection
    graph_projection_cypher = f"""
    CALL gds.graph.project(
      'QWS_it{myLabel}',
      'QWS_it{myLabel}',
      {{
        RELATIONSHIP: {{
          orientation: 'UNDIRECTED',
          type: '*',
          properties: {{
            num_common_samples: {{
              property: 'num_common_samples'
            }}
          }}
        }}
      }}
    )
    """

    # Modularity optimization community detection
    modopt_cypher = f"""
    CALL gds.modularityOptimization.stream('QWS_it{myLabel}', {{ relationshipWeightProperty: 'num_common_samples' }})
    YIELD nodeId, communityId
    RETURN nodeId, communityId,
    gds.util.asNode(nodeId).img_name AS img_name, 
    gds.util.asNode(nodeId).tissue_region AS tissue_region,
    gds.util.asNode(nodeId).event AS event,
    gds.util.asNode(nodeId).calc as calc,
    gds.util.asNode(nodeId).leaf_status AS leaf_status,
    gds.util.asNode(nodeId).median_signal AS median_signal,
    gds.util.asNode(nodeId).molecular_data AS molecular_data,
    gds.util.asNode(nodeId).no_of_samples AS no_of_samples,
    gds.util.asNode(nodeId).quartile AS quartile,
    gds.util.asNode(nodeId).samples AS samples,
    gds.util.asNode(nodeId).subsample AS subsample
    ORDER BY communityId DESC;
    """

    csv_file_path = f"/Users/m254284/Desktop/Neo4j_Projects/Img_events/ModOpt_bycalc/Iterations/QWS/modopt_calc_it{myLabel}.csv"

    # Open the CSV file for writing
    with open(csv_file_path, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        # Write header row
        csv_writer.writerow(['node_id','community','img_name', 'tissue_region','event','calc','leaf_status',
                             'median_signal', 'molecular_data','no_of_samples', 'quartile', 'samples', 'subsample'])

        with GraphDatabase.driver(uri, auth=(username, password)) as driver:
            with driver.session() as session:
                # Add label to selected nodes
                session.run(add_label_cypher)

                # Create graph projection
                session.run(graph_projection_cypher)

                # Run modularity optimization community detection
                result = session.run(modopt_cypher)

                # Write the result to the CSV file
                for record in result:
                    node_id = record['nodeId']
                    community = record['communityId']
                    img_name = record['img_name']
                    tissue_region = record['tissue_region']
                    event = record['event']
                    calc = record['calc']
                    leaf_status = record['leaf_status']
                    median_signal = record['median_signal']
                    molecular_data = record['molecular_data']
                    no_of_samples = record['no_of_samples']
                    quartile = record['quartile']
                    samples = record['samples']
                    subsample = record['subsample']

                    csv_writer.writerow([node_id, community, img_name, tissue_region, event, calc, leaf_status,
                                     median_signal, molecular_data, no_of_samples, quartile, samples, subsample])

    # Close the Neo4j driver
    driver.close()

    # Print a message indicating completion for each label
    print(f"Finished all communities in QWS_it{myLabel}.")


# ---------- UCLA

# Iterate through each label in BSW
for myLabel in labels:
    # Add label to selected nodes
    add_label_cypher = f"""
    MATCH (n)
    WHERE n.no_of_samples <> 1 AND n.no_of_samples <> 103 AND (n.calc = 'UCLA' OR n.calc = 'standard')
    WITH n.img_name AS imgName, COLLECT(n) AS nodes
    WITH imgName, nodes, nodes[TOINTEGER(RAND() * SIZE(nodes))] AS randomNode
    SET randomNode:UCLA_it{myLabel}
    """

    # Create graph projection
    graph_projection_cypher = f"""
    CALL gds.graph.project(
      'UCLA_it{myLabel}',
      'UCLA_it{myLabel}',
      {{
        RELATIONSHIP: {{
          orientation: 'UNDIRECTED',
          type: '*',
          properties: {{
            num_common_samples: {{
              property: 'num_common_samples'
            }}
          }}
        }}
      }}
    )
    """

    # Modularity optimization community detection
    modopt_cypher = f"""
    CALL gds.modularityOptimization.stream('UCLA_it{myLabel}', {{ relationshipWeightProperty: 'num_common_samples' }})
    YIELD nodeId, communityId
    RETURN nodeId, communityId,
    gds.util.asNode(nodeId).img_name AS img_name, 
    gds.util.asNode(nodeId).tissue_region AS tissue_region,
    gds.util.asNode(nodeId).event AS event,
    gds.util.asNode(nodeId).calc as calc,
    gds.util.asNode(nodeId).leaf_status AS leaf_status,
    gds.util.asNode(nodeId).median_signal AS median_signal,
    gds.util.asNode(nodeId).molecular_data AS molecular_data,
    gds.util.asNode(nodeId).no_of_samples AS no_of_samples,
    gds.util.asNode(nodeId).quartile AS quartile,
    gds.util.asNode(nodeId).samples AS samples,
    gds.util.asNode(nodeId).subsample AS subsample
    ORDER BY communityId DESC;
    """

    csv_file_path = f"/Users/m254284/Desktop/Neo4j_Projects/Img_events/ModOpt_bycalc/Iterations/UCLA/modopt_calc_it{myLabel}.csv"

    # Open the CSV file for writing
    with open(csv_file_path, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        # Write header row
        csv_writer.writerow(['node_id','community','img_name', 'tissue_region','event','calc','leaf_status',
                             'median_signal', 'molecular_data','no_of_samples', 'quartile', 'samples', 'subsample'])

        with GraphDatabase.driver(uri, auth=(username, password)) as driver:
            with driver.session() as session:
                # Add label to selected nodes
                session.run(add_label_cypher)

                # Create graph projection
                session.run(graph_projection_cypher)

                # Run modularity optimization community detection
                result = session.run(modopt_cypher)

                # Write the result to the CSV file
                for record in result:
                    node_id = record['nodeId']
                    community = record['communityId']
                    img_name = record['img_name']
                    tissue_region = record['tissue_region']
                    event = record['event']
                    calc = record['calc']
                    leaf_status = record['leaf_status']
                    median_signal = record['median_signal']
                    molecular_data = record['molecular_data']
                    no_of_samples = record['no_of_samples']
                    quartile = record['quartile']
                    samples = record['samples']
                    subsample = record['subsample']

                    csv_writer.writerow([node_id, community, img_name, tissue_region, event, calc, leaf_status,
                                     median_signal, molecular_data, no_of_samples, quartile, samples, subsample])

    # Close the Neo4j driver
    driver.close()

    # Print a message indicating completion for each label
    print(f"Finished all communities in UCLA_it{myLabel}.")




# nc -----------

# Iterate through each label in BSW
for myLabel in labels:
    # Add label to selected nodes
    add_label_cypher = f"""
    MATCH (n)
    WHERE n.no_of_samples <> 1 AND n.no_of_samples <> 103 AND (n.calc = 'nc' OR n.calc = 'none' OR n.calc = 'standard')
    WITH n.img_name AS imgName, COLLECT(n) AS nodes
    WITH imgName, nodes, nodes[TOINTEGER(RAND() * SIZE(nodes))] AS randomNode
    SET randomNode:nc_it{myLabel}
    """

    # Create graph projection
    graph_projection_cypher = f"""
    CALL gds.graph.project(
      'nc_it{myLabel}',
      'nc_it{myLabel}',
      {{
        RELATIONSHIP: {{
          orientation: 'UNDIRECTED',
          type: '*',
          properties: {{
            num_common_samples: {{
              property: 'num_common_samples'
            }}
          }}
        }}
      }}
    )
    """

    # Modularity optimization community detection
    modopt_cypher = f"""
    CALL gds.modularityOptimization.stream('nc_it{myLabel}', {{ relationshipWeightProperty: 'num_common_samples' }})
    YIELD nodeId, communityId
    RETURN nodeId, communityId,
    gds.util.asNode(nodeId).img_name AS img_name, 
    gds.util.asNode(nodeId).tissue_region AS tissue_region,
    gds.util.asNode(nodeId).event AS event,
    gds.util.asNode(nodeId).calc as calc,
    gds.util.asNode(nodeId).leaf_status AS leaf_status,
    gds.util.asNode(nodeId).median_signal AS median_signal,
    gds.util.asNode(nodeId).molecular_data AS molecular_data,
    gds.util.asNode(nodeId).no_of_samples AS no_of_samples,
    gds.util.asNode(nodeId).quartile AS quartile,
    gds.util.asNode(nodeId).samples AS samples,
    gds.util.asNode(nodeId).subsample AS subsample
    ORDER BY communityId DESC;
    """

    csv_file_path = f"/Users/m254284/Desktop/Neo4j_Projects/Img_events/ModOpt_bycalc/Iterations/nc/modopt_calc_it{myLabel}.csv"

    # Open the CSV file for writing
    with open(csv_file_path, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        # Write header row
        csv_writer.writerow(['node_id','community','img_name', 'tissue_region','event','calc','leaf_status',
                             'median_signal', 'molecular_data','no_of_samples', 'quartile', 'samples', 'subsample'])

        with GraphDatabase.driver(uri, auth=(username, password)) as driver:
            with driver.session() as session:
                # Add label to selected nodes
                session.run(add_label_cypher)

                # Create graph projection
                session.run(graph_projection_cypher)

                # Run modularity optimization community detection
                result = session.run(modopt_cypher)

                # Write the result to the CSV file
                for record in result:
                    node_id = record['nodeId']
                    community = record['communityId']
                    img_name = record['img_name']
                    tissue_region = record['tissue_region']
                    event = record['event']
                    calc = record['calc']
                    leaf_status = record['leaf_status']
                    median_signal = record['median_signal']
                    molecular_data = record['molecular_data']
                    no_of_samples = record['no_of_samples']
                    quartile = record['quartile']
                    samples = record['samples']
                    subsample = record['subsample']

                    csv_writer.writerow([node_id, community, img_name, tissue_region, event, calc, leaf_status,
                                     median_signal, molecular_data, no_of_samples, quartile, samples, subsample])

    # Close the Neo4j driver
    driver.close()

    # Print a message indicating completion for each label
    print(f"Finished all communities in nc_it{myLabel}.")



