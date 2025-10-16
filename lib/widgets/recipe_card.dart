import 'package:flutter/material.dart';

/// A reusable widget that renders information about a recipe.
///
/// The previous implementation in the project relied heavily on deeply nested
/// list literals and conditional spreads, which made it easy to misplace
/// brackets. This rewritten version keeps the structure intentionally simple so
/// the layout is easier to maintain and the code is less error-prone.
class RecipeCard extends StatelessWidget {
  const RecipeCard({
    super.key,
    required this.title,
    this.subtitle,
    this.imageUrl,
    this.duration,
    this.servings,
    this.ingredients = const <String>[],
    this.ingredientDetails = const <RecipeIngredient>[],
    this.instructions = const <String>[],
    this.onTap,
  });

  /// The recipe title.
  final String title;

  /// Optional subtitle/description.
  final String? subtitle;

  /// Optional URL to an image that represents the recipe.
  final String? imageUrl;

  /// Total time required to prepare the recipe.
  final Duration? duration;

  /// Number of servings the recipe yields.
  final int? servings;

  /// Ingredients required for the recipe, represented as free-form strings.
  ///
  /// This is kept for backwards compatibility. When [ingredientDetails] is
  /// provided the structured values take precedence.
  final List<String> ingredients;

  /// Structured ingredients with dedicated fields for quantity and name.
  final List<RecipeIngredient> ingredientDetails;

  /// Step-by-step instructions for preparing the recipe.
  final List<String> instructions;

  /// Callback triggered when the card is tapped.
  final VoidCallback? onTap;

  @override
  Widget build(BuildContext context) {
    final theme = Theme.of(context);
    final textTheme = theme.textTheme;

    final structuredIngredients = _resolveIngredients();

    return Card(
      clipBehavior: Clip.antiAlias,
      child: InkWell(
        onTap: onTap,
        child: Column(
          crossAxisAlignment: CrossAxisAlignment.stretch,
          children: [
            if (imageUrl != null && imageUrl!.isNotEmpty)
              _RecipeImage(imageUrl: imageUrl!),
            Padding(
              padding: const EdgeInsets.all(16),
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text(title, style: textTheme.titleLarge),
                  if (subtitle != null && subtitle!.isNotEmpty) ...[
                    const SizedBox(height: 8),
                    Text(
                      subtitle!,
                      style: textTheme.bodyMedium?.copyWith(
                        color: theme.colorScheme.onSurfaceVariant,
                      ),
                    ),
                  ],
                  if (duration != null || servings != null) ...[
                    const SizedBox(height: 12),
                    Wrap(
                      spacing: 12,
                      runSpacing: 4,
                      children: [
                        if (duration != null)
                          _InfoChip(
                            icon: Icons.schedule,
                            label: _formatDuration(duration!),
                          ),
                        if (servings != null)
                          _InfoChip(
                            icon: Icons.restaurant,
                            label: '$servings servings',
                          ),
                      ],
                    ),
                  ],
                  if (structuredIngredients.isNotEmpty) ...[
                    const SizedBox(height: 16),
                    Text('Ingredients', style: textTheme.titleMedium),
                    const SizedBox(height: 8),
                    _IngredientList(ingredients: structuredIngredients),
                  ] else if (ingredients.isNotEmpty) ...[
                    const SizedBox(height: 16),
                    Text('Ingredients', style: textTheme.titleMedium),
                    const SizedBox(height: 8),
                    ...ingredients.map(
                      (ingredient) => Padding(
                        padding: const EdgeInsets.only(bottom: 4),
                        child: Row(
                          crossAxisAlignment: CrossAxisAlignment.start,
                          children: [
                            const Text('• '),
                            Expanded(
                              child: Text(
                                ingredient,
                                style: textTheme.bodyMedium,
                              ),
                            ),
                          ],
                        ),
                      ),
                    ),
                  ],
                  if (instructions.isNotEmpty) ...[
                    const SizedBox(height: 16),
                    Text('Instructions', style: textTheme.titleMedium),
                    const SizedBox(height: 8),
                    ...instructions.asMap().entries.map(
                      (entry) => Padding(
                        padding: const EdgeInsets.only(bottom: 8),
                        child: Row(
                          crossAxisAlignment: CrossAxisAlignment.start,
                          children: [
                            Text('${entry.key + 1}. '),
                            Expanded(
                              child: Text(
                                entry.value,
                                style: textTheme.bodyMedium,
                              ),
                            ),
                          ],
                        ),
                      ),
                    ),
                  ],
                ],
              ),
            ),
          ],
        ),
      ),
    );
  }

  List<RecipeIngredient> _resolveIngredients() {
    final fromDetails = ingredientDetails
        .where((item) => item.name.trim().isNotEmpty || item.hasQuantity)
        .toList();

    if (fromDetails.isNotEmpty) {
      return fromDetails;
    }

    return ingredients
        .map((item) => RecipeIngredient(name: item))
        .where((item) => item.name.trim().isNotEmpty)
        .toList();
  }
}

class _RecipeImage extends StatelessWidget {
  const _RecipeImage({required this.imageUrl});

  final String imageUrl;

  @override
  Widget build(BuildContext context) {
    return AspectRatio(
      aspectRatio: 16 / 9,
      child: Image.network(
        imageUrl,
        fit: BoxFit.cover,
        errorBuilder: (context, error, stackTrace) => Container(
          color: Theme.of(context).colorScheme.surfaceVariant,
          alignment: Alignment.center,
          child: const Icon(Icons.image_not_supported),
        ),
      ),
    );
  }
}

class _InfoChip extends StatelessWidget {
  const _InfoChip({
    required this.icon,
    required this.label,
  });

  final IconData icon;
  final String label;

  @override
  Widget build(BuildContext context) {
    final theme = Theme.of(context);
    return Container(
      padding: const EdgeInsets.symmetric(horizontal: 10, vertical: 6),
      decoration: BoxDecoration(
        color: theme.colorScheme.surfaceVariant,
        borderRadius: BorderRadius.circular(20),
      ),
      child: Row(
        mainAxisSize: MainAxisSize.min,
        children: [
          Icon(icon, size: 16, color: theme.colorScheme.onSurfaceVariant),
          const SizedBox(width: 6),
          Text(
            label,
            style: theme.textTheme.labelMedium?.copyWith(
              color: theme.colorScheme.onSurfaceVariant,
            ),
          ),
        ],
      ),
    );
  }
}

class _IngredientList extends StatelessWidget {
  const _IngredientList({required this.ingredients});

  final List<RecipeIngredient> ingredients;

  static const double _quantityColumnWidth = 96;

  @override
  Widget build(BuildContext context) {
    final theme = Theme.of(context);
    final textTheme = theme.textTheme;

    return Container(
      width: double.infinity,
      decoration: BoxDecoration(
        color: theme.colorScheme.surfaceVariant.withOpacity(0.4),
        borderRadius: BorderRadius.circular(12),
      ),
      padding: const EdgeInsets.symmetric(horizontal: 12, vertical: 12),
      child: Column(
        children: [
          for (final ingredient in ingredients)
            Padding(
              padding: const EdgeInsets.symmetric(vertical: 6),
              child: _IngredientRow(
                ingredient: ingredient,
                quantityWidth: _quantityColumnWidth,
                textTheme: textTheme,
                theme: theme,
              ),
            ),
        ],
      ),
    );
  }
}

class _IngredientRow extends StatelessWidget {
  const _IngredientRow({
    required this.ingredient,
    required this.quantityWidth,
    required this.textTheme,
    required this.theme,
  });

  final RecipeIngredient ingredient;
  final double quantityWidth;
  final TextTheme textTheme;
  final ThemeData theme;

  @override
  Widget build(BuildContext context) {
    final hasQuantity = ingredient.hasQuantity;

    return Row(
      crossAxisAlignment: CrossAxisAlignment.start,
      children: [
        Container(
          margin: const EdgeInsets.only(top: 6),
          width: 6,
          height: 6,
          decoration: BoxDecoration(
            color: theme.colorScheme.primary,
            shape: BoxShape.circle,
          ),
        ),
        const SizedBox(width: 12),
        SizedBox(
          width: quantityWidth,
          child: hasQuantity
              ? Text(
                  ingredient.quantity!,
                  style: textTheme.bodyMedium?.copyWith(
                    fontWeight: FontWeight.w600,
                    color: theme.colorScheme.onSurface,
                  ),
                  textAlign: TextAlign.right,
                )
              : const SizedBox.shrink(),
        ),
        const SizedBox(width: 12),
        Expanded(
          child: Text(
            ingredient.name,
            style: textTheme.bodyMedium,
          ),
        ),
      ],
    );
  }
}

/// Representation of a single ingredient with optional quantity information.
class RecipeIngredient {
  const RecipeIngredient({
    required this.name,
    this.quantity,
  });

  /// Creates an ingredient from a JSON-like map.
  factory RecipeIngredient.fromMap(Map<String, dynamic> map) {
    return RecipeIngredient(
      name: (map['name'] ?? '').toString(),
      quantity: map['quantity']?.toString(),
    );
  }

  final String name;
  final String? quantity;

  bool get hasQuantity => quantity != null && quantity!.trim().isNotEmpty;
}

String _formatDuration(Duration duration) {
  final hours = duration.inHours;
  final minutes = duration.inMinutes.remainder(60);
  if (hours > 0 && minutes > 0) {
    return '$hours h $minutes min';
  }
  if (hours > 0) {
    return '$hours h';
  }
  return '$minutes min';
}
